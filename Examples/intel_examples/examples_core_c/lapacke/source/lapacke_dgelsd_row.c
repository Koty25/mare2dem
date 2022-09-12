/*******************************************************************************
* Copyright 2010-2020 Intel Corporation.
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
   LAPACKE_dgelsd Example.
   =======================

   Program computes the minimum norm-solution to a real linear least squares
   problem using the singular value decomposition of A,
   where A is the coefficient matrix:

     0.12  -8.19   7.69  -2.26  -4.71
    -6.91   2.22  -5.12  -9.08   9.96
    -3.33  -8.94  -6.72  -4.40  -9.98
     3.97   3.33  -2.74  -7.92  -3.20

   and B is the right-hand side matrix:

     7.30   0.47  -6.28
     1.33   6.58  -3.42
     2.68  -1.71   3.46
    -9.62  -0.79   0.41

   Description.
   ============

   The routine computes the minimum-norm solution to a real linear least
   squares problem: minimize ||b - A*x|| using the singular value
   decomposition (SVD) of A. A is an m-by-n matrix which may be rank-deficient.

   Several right hand side vectors b and solution vectors x can be handled
   in a single call; they are stored as the columns of the m-by-nrhs right
   hand side matrix B and the n-by-nrhs solution matrix X.

   The effective rank of A is determined by treating as zero those singular
   values which are less than rcond times the largest singular value.

   Example Program Results.
   ========================

 LAPACKE_dgelsd (row-major, high-level) Example Program Results

 Minimum norm solution
  -0.69  -0.24   0.06
  -0.80  -0.08   0.21
   0.38   0.12  -0.65
   0.29  -0.24   0.42
   0.29   0.35  -0.30

 Effective rank =      4

 Singular values
  18.66  15.99  10.01   8.51
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
#define M 4
#define N 5
#define NRHS 3
#define LDA N
#define LDB NRHS

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond = -1.0;
	/* Local arrays */
	double s[M];
	double a[LDA*M] = {
	    0.12, -8.19, 7.69, -2.26, -4.71,
	   -6.91,  2.22, -5.12, -9.08,  9.96,
	   -3.33, -8.94, -6.72, -4.40, -9.98,
	    3.97,  3.33, -2.74, -7.92, -3.20
	};
	double b[LDB*N] = {
	    7.30,  0.47, -6.28,
	    1.33,  6.58, -3.42,
	    2.68, -1.71, 3.46,
	   -9.62, -0.79, 0.41,
	    0.00,  0.00, 0.00
	};
	/* Executable statements */
	printf( "LAPACKE_dgelsd (row-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */
	info = LAPACKE_dgelsd( LAPACK_ROW_MAJOR, m, n, nrhs, a, lda, b, ldb,
			s, rcond, &rank );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm computing SVD failed to converge;\n" );
		printf( "the least squares solution could not be computed.\n" );
		exit( 1 );
	}
	/* Print minimum norm solution */
	print_matrix( "Minimum norm solution", n, nrhs, b, ldb );
	/* Print effective rank */
	printf( "\n Effective rank = %6i\n", rank );
	/* Print singular values */
	print_matrix( "Singular values", 1, m, s, 1 );
	exit( 0 );
} /* End of LAPACKE_dgelsd Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
		printf( "\n" );
	}
}
