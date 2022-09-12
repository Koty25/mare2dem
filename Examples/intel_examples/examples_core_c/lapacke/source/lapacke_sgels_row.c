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
   LAPACKE_sgels Example.
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

 LAPACKE_sgels (row-major, high-level) Example Program Results

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
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda );
extern void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda );

/* Parameters */
#define M 6
#define N 4
#define NRHS 2
#define LDA N
#define LDB NRHS

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
	/* Local arrays */
	float a[LDA*M] = {
	    1.44f, -7.84f, -4.39f,  4.53f,
	   -9.96f, -0.28f, -3.24f,  3.83f,
	   -7.55f, 3.24f, 6.27f, -6.64f,
	    8.34f, 8.09f, 5.28f,  2.06f,
	    7.08f, 2.52f, 0.74f, -2.47f,
	   -5.45f, -5.70f, -1.19f,  4.70f
	};
	float b[LDB*M] = {
	    8.58f, 9.35f,
	    8.26f, -4.43f,
	    8.48f, -0.70f,
	   -5.28f, -0.26f,
	    5.72f, -7.36f,
	    8.93f, -2.52f
	};
	/* Executable statements */
	printf( "LAPACKE_sgels (row-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */
	info = LAPACKE_sgels( LAPACK_ROW_MAJOR, 'N', m, n, nrhs, a, lda,
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
			&b[n*ldb], ldb );
	/* Print details of QR factorization */
	print_matrix( "Details of QR factorization", m, n, a, lda );
	exit( 0 );
} /* End of LAPACKE_sgels Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing norms of matrix columns */
void print_vector_norm( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda ) {
	MKL_INT i, j;
	float norm;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) {
		norm = 0.0;
		for( i = 0; i < m; i++ ) norm += a[i*lda+j] * a[i*lda+j];
		printf( " %6.2f", norm );
	}
	printf( "\n" );
}
