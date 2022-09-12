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
   CHEEVD Example.
   ==============

   Program computes all eigenvalues and eigenvectors of a complex Hermitian
   matrix A using divide and conquer algorithm, where A is:

   (  3.40,  0.00) ( -2.36, -1.93) ( -4.68,  9.55) (  5.37, -1.23)
   ( -2.36,  1.93) (  6.94,  0.00) (  8.13, -1.47) (  2.07, -5.78)
   ( -4.68, -9.55) (  8.13,  1.47) ( -2.14,  0.00) (  4.68,  7.44)
   (  5.37,  1.23) (  2.07,  5.78) (  4.68, -7.44) ( -7.42,  0.00)

   Description.
   ============

   The routine computes all eigenvalues and, optionally, eigenvectors of an
   n-by-n complex Hermitian matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.
   If the eigenvectors are requested, then this routine uses a divide and
   conquer algorithm to compute eigenvalues and eigenvectors.

   Example Program Results.
   ========================

 CHEEVD Example Program Results

 Eigenvalues
 -21.97  -0.05   6.46  16.34

 Eigenvectors (stored columnwise)
 (  0.41,  0.00) ( -0.34,  0.00) ( -0.69,  0.00) (  0.49,  0.00)
 (  0.02, -0.30) (  0.32, -0.21) ( -0.57, -0.22) ( -0.59, -0.21)
 (  0.18,  0.57) ( -0.42, -0.32) (  0.06,  0.16) ( -0.35, -0.47)
 ( -0.62, -0.09) ( -0.58,  0.35) ( -0.15, -0.31) ( -0.10, -0.12)
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
#define N 4
#define LDA N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, lda = LDA, info, lwork, lrwork, liwork;
	MKL_INT iwkopt;
	MKL_INT* iwork;
	float rwkopt;
	float* rwork;
	fcomplex wkopt;
	fcomplex* work;
	/* Local arrays */
	float w[N];
	fcomplex a[LDA*N] = {
	   { 3.40f,  0.00f}, {-2.36f,  1.93f}, {-4.68f, -9.55f}, { 5.37f,  1.23f},
	   { 0.00f,  0.00f}, { 6.94f,  0.00f}, { 8.13f,  1.47f}, { 2.07f,  5.78f},
	   { 0.00f,  0.00f}, { 0.00f,  0.00f}, {-2.14f,  0.00f}, { 4.68f, -7.44f},
	   { 0.00f,  0.00f}, { 0.00f,  0.00f}, { 0.00f,  0.00f}, {-7.42f,  0.00f}
	};
	/* Executable statements */
	printf( " CHEEVD Example Program Results\n" );
	/* Query and allocate the optimal workspace */
	lwork = -1;
	lrwork = -1;
	liwork = -1;
	cheevd( "Vectors", "Lower", &n, a, &lda, w, &wkopt, &lwork, &rwkopt,
			&lrwork, &iwkopt, &liwork, &info );
	lwork = (MKL_INT)wkopt.re;
	work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
	lrwork = (MKL_INT)rwkopt;
	rwork = (float*)malloc( lrwork*sizeof(float) );
	liwork = iwkopt;
	iwork = (MKL_INT*)malloc( liwork*sizeof(MKL_INT) );
	/* Solve eigenproblem */
	cheevd( "Vectors", "Lower", &n, a, &lda, w, work, &lwork, rwork,
			&lrwork, iwork, &liwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	/* Print eigenvalues */
	print_rmatrix( "Eigenvalues", 1, n, w, 1 );
	/* Print eigenvectors */
	print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
	/* Free workspace */
	free( (void*)iwork );
	free( (void*)rwork );
	free( (void*)work );
	exit( 0 );
} /* End of CHEEVD Example */

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
