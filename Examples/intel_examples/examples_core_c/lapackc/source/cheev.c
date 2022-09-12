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
   CHEEV Example.
   ==============

   Program computes all eigenvalues and eigenvectors of a complex Hermitian
   matrix A:

   (  9.14,  0.00) ( -4.37, -9.22) ( -1.98, -1.72) ( -8.96, -9.50)
   ( -4.37,  9.22) ( -3.35,  0.00) (  2.25, -9.51) (  2.57,  2.40)
   ( -1.98,  1.72) (  2.25,  9.51) ( -4.82,  0.00) ( -3.24,  2.04)
   ( -8.96,  9.50) (  2.57, -2.40) ( -3.24, -2.04) (  8.44,  0.00)

   Description.
   ============

   The routine computes all eigenvalues and, optionally, eigenvectors of an
   n-by-n complex Hermitian matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.

   Example Program Results.
   ========================

 CHEEV Example Program Results

 Eigenvalues
 -16.00  -6.76   6.67  25.51

 Eigenvectors (stored columnwise)
 (  0.34,  0.00) ( -0.55,  0.00) (  0.31,  0.00) ( -0.70,  0.00)
 (  0.44, -0.54) (  0.26,  0.18) (  0.45,  0.29) (  0.22, -0.28)
 ( -0.48, -0.37) ( -0.52, -0.02) ( -0.05,  0.57) (  0.15,  0.08)
 (  0.10, -0.12) ( -0.50,  0.28) ( -0.23, -0.48) (  0.34, -0.49)
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
	MKL_INT n = N, lda = LDA, info, lwork;
	fcomplex wkopt;
	fcomplex* work;
	/* Local arrays */
	/* rwork dimension should be at least max(1,3*n-2) */
	float w[N], rwork[3*N-2];
	fcomplex a[LDA*N] = {
	   { 9.14f,  0.00f}, {-4.37f,  9.22f}, {-1.98f,  1.72f}, {-8.96f,  9.50f},
	   { 0.00f,  0.00f}, {-3.35f,  0.00f}, { 2.25f,  9.51f}, { 2.57f, -2.40f},
	   { 0.00f,  0.00f}, { 0.00f,  0.00f}, {-4.82f,  0.00f}, {-3.24f, -2.04f},
	   { 0.00f,  0.00f}, { 0.00f,  0.00f}, { 0.00f,  0.00f}, { 8.44f,  0.00f}
	};
	/* Executable statements */
	printf( " CHEEV Example Program Results\n" );
	/* Query and allocate the optimal workspace */
	lwork = -1;
	cheev( "Vectors", "Lower", &n, a, &lda, w, &wkopt, &lwork, rwork, &info );
	lwork = (MKL_INT)wkopt.re;
	work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
	/* Solve eigenproblem */
	cheev( "Vectors", "Lower", &n, a, &lda, w, work, &lwork, rwork, &info );
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
	free( (void*)work );
	exit( 0 );
} /* End of CHEEV Example */

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
