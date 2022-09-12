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
   LAPACKE_zgelsd Example.
   =======================

   Program computes the minimum norm-solution to a complex linear least squares
   problem using the singular value decomposition of A,
   where A is the coefficient matrix:

   (  4.55, -0.32) ( -4.36, -4.76) (  3.99, -6.84) (  8.03, -6.47)
   (  8.87, -3.11) (  0.02,  8.43) (  5.43, -9.30) (  2.28,  8.94)
   ( -0.74,  1.16) (  3.80, -6.12) ( -7.24,  0.72) (  2.21,  9.52)

   and B is the right-hand side matrix:

   ( -8.25,  7.98) (  2.91, -8.81)
   ( -5.04,  3.33) (  6.19,  0.19)
   (  7.98, -4.38) ( -5.96,  7.18)

   Description.
   ============

   The routine computes the minimum-norm solution to a complex linear least
   squares problem: minimize ||b - A*x|| using the singular value
   decomposition (SVD) of A. A is an m-by-n matrix which may be rank-deficient.

   Several right hand side vectors b and solution vectors x can be handled
   in a single call; they are stored as the columns of the m-by-nrhs right
   hand side matrix B and the n-by-nrhs solution matrix X.

   The effective rank of A is determined by treating as zero those singular
   values which are less than rcond times the largest singular value.

   Example Program Results.
   ========================

 LAPACKE_zgelsd (column-major, high-level) Example Program Results

 Minimum norm solution
 ( -0.08,  0.09) (  0.04,  0.16)
 ( -0.17,  0.10) (  0.17, -0.47)
 ( -0.92, -0.01) (  0.71, -0.41)
 ( -0.47, -0.26) (  0.69,  0.02)

 Effective rank =      3

 Singular values
  20.01  18.21   7.88
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda );
extern void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
#define M 3
#define N 4
#define NRHS 2
#define LDA M
#define LDB N

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, rank;
	/* Negative rcond means using default (machine precision) value */
	double rcond = -1.0;
	/* Local arrays */
	double s[M];
	MKL_Complex16 a[LDA*N] = {
	   { 4.55, -0.32}, { 8.87, -3.11}, {-0.74,  1.16},
	   {-4.36, -4.76}, { 0.02,  8.43}, { 3.80, -6.12},
	   { 3.99, -6.84}, { 5.43, -9.30}, {-7.24,  0.72},
	   { 8.03, -6.47}, { 2.28,  8.94}, { 2.21,  9.52}
	};
	MKL_Complex16 b[LDB*NRHS] = {
	   {-8.25,  7.98}, {-5.04,  3.33}, { 7.98, -4.38}, { 0.00,  0.00},
	   { 2.91, -8.81}, { 6.19,  0.19}, {-5.96,  7.18}, { 0.00,  0.00}
	};
	/* Executable statements */
	printf( "LAPACKE_zgelsd (column-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */
	info = LAPACKE_zgelsd( LAPACK_COL_MAJOR, m, n, nrhs, a, lda, b, ldb,
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
	print_rmatrix( "Singular values", 1, m, s, 1 );
	exit( 0 );
} /* End of LAPACKE_zgelsd Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%6.2f,%6.2f)", a[i+j*lda].real, a[i+j*lda].imag );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}
