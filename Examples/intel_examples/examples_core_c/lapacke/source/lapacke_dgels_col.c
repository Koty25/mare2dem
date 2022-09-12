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
   LAPACKE_dgels Example.
   ======================

   Program computes the least squares solution to the overdetermined linear
   system A*X = B with full rank matrix A using QR factorization,
   where A is the coefficient matrix:

     1.44  -7.84  -4.39   4.53
    -9.96  -0.28  -3.24   3.83
    -7.55   3.24   6.27  -6.64
     8.34   8.09   5.28   2.06
     7.08   2.52   0.74  -2.47
    -5.45  -5.70  -1.19   4.70

   and B is the right-hand side matrix:

     8.58   9.35
     8.26  -4.43
     8.48  -0.70
    -5.28  -0.26
     5.72  -7.36
     8.93  -2.52

   Description.
   ============

   The routine solves overdetermined or underdetermined real linear systems
   involving an m-by-n matrix A, or its transpose, using a QR or LQ
   factorization of A. It is assumed that A has full rank.

   Several right hand side vectors b and solution vectors x can be handled
   in a single call; they are stored as the columns of the m-by-nrhs right
   hand side matrix B and the n-by-nrhs solution matrix X.

   Example Program Results.
   ========================

 LAPACKE_dgels (column-major, high-level) Example Program Results

 Solution
  -0.45   0.25
  -0.85  -0.90
   0.71   0.63
   0.13   0.14

 Residual sum of squares for the solution
 195.36 107.06

 Details of QR factorization
 -17.54  -4.76  -1.96   0.42
  -0.52  12.40   7.88  -5.84
  -0.40  -0.14  -5.75   4.11
   0.44  -0.66  -0.20  -7.78
   0.37  -0.26  -0.17  -0.15
  -0.29   0.46   0.41   0.24
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );
extern void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
#define M 6
#define N 4
#define NRHS 2
#define LDA M
#define LDB M

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
	/* Local arrays */
	double a[LDA*N] = {
	    1.44, -9.96, -7.55,  8.34,  7.08, -5.45,
	   -7.84, -0.28,  3.24,  8.09,  2.52, -5.70,
	   -4.39, -3.24,  6.27,  5.28,  0.74, -1.19,
	    4.53,  3.83, -6.64,  2.06, -2.47,  4.70
	};
	double b[LDB*NRHS] = {
	    8.58,  8.26,  8.48, -5.28,  5.72,  8.93,
	    9.35, -4.43, -0.70, -0.26, -7.36, -2.52
	};
	/* Executable statements */
	printf( "LAPACKE_dgels (column-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */
	info = LAPACKE_dgels( LAPACK_COL_MAJOR, 'N', m, n, nrhs, a, lda,
			b, ldb );
	/* Check for the full rank */
	if( info > 0 ) {
		printf( "The diagonal element %i of the triangular factor ", info );
		printf( "of A is zero, so that A does not have full rank;\n" );
		printf( "the least squares solution could not be computed.\n" );
		exit( 1 );
	}
	/* Print least squares solution */
	print_matrix( "Least squares solution", n, nrhs, b, ldb );
	/* Print residual sum of squares for the solution */
	print_vector_norm( "Residual sum of squares for the solution", m-n, nrhs,
			&b[n], ldb );
	/* Print details of QR factorization */
	print_matrix( "Details of QR factorization", m, n, a, lda );
	exit( 0 );
} /* End of LAPACKE_dgels Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing norms of matrix columns */
void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	double norm;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) {
		norm = 0.0;
		for( i = 0; i < m; i++ ) norm += a[i+j*lda] * a[i+j*lda];
		printf( " %6.2f", norm );
	}
	printf( "\n" );
}
