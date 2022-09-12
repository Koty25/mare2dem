/*******************************************************************************
* Copyright 2009-2020 Intel Corporation.
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
   CGELSD Example.
   ==============

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

 CGELSD Example Program Results

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
#include <mkl.h>

/* Complex datatype */
struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, fcomplex* a, int lda );
extern void print_rmatrix( char* desc, int m, int n, float* a, int lda );

/* Parameters */
#define M 3
#define N 4
#define NRHS 2
#define LDA M
#define LDB N

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, lwork, rank;
	/* Negative rcond means using default (machine precision) value */
	float rcond = -1.0;
	fcomplex wkopt;
	fcomplex* work;
	/* Local arrays */
	/* iwork dimension should be at least 3*min(m,n)*nlvl + 11*min(m,n),
		rwork dimension should be at least 10*min(m,n)+2*min(m,n)*smlsiz+
		+8*min(m,n)*nlvl+3*smlsiz*nrhs+(smlsiz+1)^2,
		where nlvl = max( 0, int( log_2( min(m,n)/(smlsiz+1) ) )+1 )
		and smlsiz = 25 */
	MKL_INT iwork[3*M*0+11*M];
	float s[M], rwork[10*M+2*M*25+8*M*0+3*25*NRHS+26*26];
	fcomplex a[LDA*N] = {
	   { 4.55f, -0.32f}, { 8.87f, -3.11f}, {-0.74f,  1.16f},
	   {-4.36f, -4.76f}, { 0.02f,  8.43f}, { 3.80f, -6.12f},
	   { 3.99f, -6.84f}, { 5.43f, -9.30f}, {-7.24f,  0.72f},
	   { 8.03f, -6.47f}, { 2.28f,  8.94f}, { 2.21f,  9.52f}
	};
	fcomplex b[LDB*NRHS] = {
	   {-8.25f,  7.98f}, {-5.04f,  3.33f}, { 7.98f, -4.38f}, { 0.00f,  0.00f},
	   { 2.91f, -8.81f}, { 6.19f,  0.19f}, {-5.96f,  7.18f}, { 0.00f,  0.00f}
	};
	/* Executable statements */
	printf( " CGELSD Example Program Results\n" );
	/* Query and allocate the optimal workspace */
	lwork = -1;
	cgelsd( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, &wkopt, &lwork,
			rwork, iwork, &info );
	lwork = (MKL_INT)wkopt.re;
	work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
	/* Solve the equations A*X = B */
	cgelsd( &m, &n, &nrhs, a, &lda, b, &ldb, s, &rcond, &rank, work, &lwork,
			rwork, iwork, &info );
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
	/* Free workspace */
	free( (void*)work );
	exit( 0 );
} /* End of CGELSD Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, fcomplex* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%6.2f,%6.2f)", a[i+j*lda].re, a[i+j*lda].im );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, float* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}
