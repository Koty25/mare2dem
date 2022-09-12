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
   CGEEV Example.
   ==============

   Program computes the eigenvalues and left and right eigenvectors of a general
   rectangular matrix A:

   ( -3.84,  2.25) ( -8.94, -4.75) (  8.95, -6.53) ( -9.87,  4.82)
   ( -0.66,  0.83) ( -4.40, -3.82) ( -3.50, -4.26) ( -3.15,  7.36)
   ( -3.99, -4.73) ( -5.88, -6.60) ( -3.36, -0.40) ( -0.75,  5.23)
   (  7.74,  4.18) (  3.66, -7.53) (  2.58,  3.60) (  4.59,  5.41)

   Description.
   ============

   The routine computes for an n-by-n complex nonsymmetric matrix A, the
   eigenvalues and, optionally, the left and/or right eigenvectors. The right
   eigenvector v(j) of A satisfies

   A*v(j)= lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The left eigenvector u(j) of A satisfies

   u(j)H*A = lambda(j)*u(j)H

   where u(j)H denotes the conjugate transpose of u(j). The computed
   eigenvectors are normalized to have Euclidean norm equal to 1 and
   largest component real.

   Example Program Results.
   ========================

 CGEEV Example Program Results

 Eigenvalues
 ( -9.43,-12.98) ( -3.44, 12.69) (  0.11, -3.40) (  5.76,  7.13)

 Left eigenvectors
 (  0.24, -0.18) (  0.61,  0.00) ( -0.18, -0.33) (  0.28,  0.09)
 (  0.79,  0.00) ( -0.05, -0.27) (  0.82,  0.00) ( -0.55,  0.16)
 (  0.22, -0.27) ( -0.21,  0.53) ( -0.37,  0.15) (  0.45,  0.09)
 ( -0.02,  0.41) (  0.40, -0.24) (  0.06,  0.12) (  0.62,  0.00)

 Right eigenvectors
 (  0.43,  0.33) (  0.83,  0.00) (  0.60,  0.00) ( -0.31,  0.03)
 (  0.51, -0.03) (  0.08, -0.25) ( -0.40, -0.20) (  0.04,  0.34)
 (  0.62,  0.00) ( -0.25,  0.28) ( -0.09, -0.48) (  0.36,  0.06)
 ( -0.23,  0.11) ( -0.10, -0.32) ( -0.43,  0.13) (  0.81,  0.00)
*/
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>

/* Complex datatype */
struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, fcomplex* a, int lda );

/* Parameters */
#define N 4
#define LDA N
#define LDVL N
#define LDVR N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, lda = LDA, ldvl = LDVL, ldvr = LDVR, info, lwork;
	fcomplex wkopt;
	fcomplex* work;
	/* Local arrays */
	/* rwork dimension should be at least 2*n */
	float rwork[2*N];
	fcomplex w[N], vl[LDVL*N], vr[LDVR*N];
	fcomplex a[LDA*N] = {
	   {-3.84f,  2.25f}, {-0.66f,  0.83f}, {-3.99f, -4.73f}, { 7.74f,  4.18f},
	   {-8.94f, -4.75f}, {-4.40f, -3.82f}, {-5.88f, -6.60f}, { 3.66f, -7.53f},
	   { 8.95f, -6.53f}, {-3.50f, -4.26f}, {-3.36f, -0.40f}, { 2.58f,  3.60f},
	   {-9.87f,  4.82f}, {-3.15f,  7.36f}, {-0.75f,  5.23f}, { 4.59f,  5.41f}
	};
	/* Executable statements */
	printf( " CGEEV Example Program Results\n" );
	/* Query and allocate the optimal workspace */
	lwork = -1;
	cgeev( "Vectors", "Vectors", &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
         &wkopt, &lwork, rwork, &info );
	lwork = (MKL_INT)wkopt.re;
	work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
	/* Solve eigenproblem */
	cgeev( "Vectors", "Vectors", &n, a, &lda, w, vl, &ldvl, vr, &ldvr,
         work, &lwork, rwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	/* Print eigenvalues */
	print_matrix( "Eigenvalues", 1, n, w, 1 );
	/* Print left eigenvectors */
	print_matrix( "Left eigenvectors", n, n, vl, ldvl );
	/* Print right eigenvectors */
	print_matrix( "Right eigenvectors", n, n, vr, ldvr );
	/* Free workspace */
	free( (void*)work );
	exit( 0 );
} /* End of CGEEV Example */

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
