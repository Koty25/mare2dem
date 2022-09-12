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
   ZHEEVR Example.
   ==============

   Program computes eigenvalues specified by a selected range of values
   and corresponding eigenvectors of a complex Hermitian matrix A using the
   Relatively Robust Representations, where A is:

   ( -2.16,  0.00) ( -0.16, -4.86) ( -7.23, -9.38) ( -0.04,  6.86)
   ( -0.16,  4.86) (  7.45,  0.00) (  4.39,  6.29) ( -8.11, -4.41)
   ( -7.23,  9.38) (  4.39, -6.29) ( -9.03,  0.00) ( -6.89, -7.66)
   ( -0.04, -6.86) ( -8.11,  4.41) ( -6.89,  7.66) (  7.76,  0.00)

   Description.
   ============

   The routine computes selected eigenvalues and, optionally, eigenvectors of
   an n-by-n complex Hermitian matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.
   Eigenvalues and eigenvectors can be selected by specifying either a range
   of values or a range of indices for the desired eigenvalues.

   Example Program Results.
   ========================

 ZHEEVR Example Program Results

 The total number of eigenvalues found: 2

 Selected eigenvalues
  -4.18   3.57

 Selected eigenvectors (stored columnwise)
 (  0.68,  0.00) (  0.38,  0.00)
 (  0.03,  0.18) (  0.54, -0.57)
 ( -0.03,  0.21) ( -0.40,  0.04)
 (  0.20,  0.64) ( -0.14, -0.26)
*/
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>

/* Complex datatype */
struct _dcomplex { double re, im; };
typedef struct _dcomplex dcomplex;

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, dcomplex* a, int lda );
extern void print_rmatrix( char* desc, int m, int n, double* a, int lda );

/* Parameters */
#define N 4
#define LDA N
#define LDZ N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, lda = LDA, ldz = LDZ, il, iu, m, info, lwork, lrwork, liwork;
	double abstol, vl, vu;
	MKL_INT iwkopt;
	MKL_INT* iwork;
	double rwkopt;
	double* rwork;
	dcomplex wkopt;
	dcomplex* work;
	/* Local arrays */
	MKL_INT isuppz[N];
	double w[N];
	dcomplex z[LDZ*N];
	dcomplex a[LDA*N] = {
	   {-2.16,  0.00}, {-0.16,  4.86}, {-7.23,  9.38}, {-0.04, -6.86},
	   { 0.00,  0.00}, { 7.45,  0.00}, { 4.39, -6.29}, {-8.11,  4.41},
	   { 0.00,  0.00}, { 0.00,  0.00}, {-9.03,  0.00}, {-6.89,  7.66},
	   { 0.00,  0.00}, { 0.00,  0.00}, { 0.00,  0.00}, { 7.76,  0.00}
	};
	/* Executable statements */
	printf( " ZHEEVR Example Program Results\n" );
	/* Negative abstol means using the default value */
	abstol = -1.0;
	/* Set VL, VU to compute eigenvalues in half-open (VL,VU] interval */
	vl = -5.0;
	vu = 5.0;
	/* Query and allocate the optimal workspace */
	lwork = -1;
	lrwork = -1;
	liwork = -1;
	zheevr( "Vectors", "Values", "Lower", &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &m, w, z, &ldz, isuppz, &wkopt, &lwork, &rwkopt, &lrwork,
			&iwkopt, &liwork, &info );
	lwork = (MKL_INT)wkopt.re;
	work = (dcomplex*)malloc( lwork*sizeof(dcomplex) );
	lrwork = (MKL_INT)rwkopt;
	rwork = (double*)malloc( lrwork*sizeof(double) );
	liwork = iwkopt;
	iwork = (MKL_INT*)malloc( liwork*sizeof(MKL_INT) );
	/* Solve eigenproblem */
	zheevr( "Vectors", "Values", "Lower", &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &m, w, z, &ldz, isuppz, work, &lwork, rwork, &lrwork,
			iwork, &liwork, &info );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	/* Print the number of eigenvalues found */
	printf( "\n The total number of eigenvalues found:%2i\n", m );
	/* Print eigenvalues */
	print_rmatrix( "Selected eigenvalues", 1, m, w, 1 );
	/* Print eigenvectors */
	print_matrix( "Selected eigenvectors (stored columnwise)", n, m, z, ldz );
	/* Free workspace */
	free( (void*)iwork );
	free( (void*)rwork );
	free( (void*)work );
	exit( 0 );
} /* End of ZHEEVR Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, dcomplex* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%6.2f,%6.2f)", a[i+j*lda].re, a[i+j*lda].im );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, int m, int n, double* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}