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
   LAPACKE_zheev Example.
   ======================

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

 LAPACKE_zheev (row-major, high-level) Example Program Results

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
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda );
extern void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
#define N 4
#define LDA N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, lda = LDA, info;
	/* Local arrays */
	double w[N];
	MKL_Complex16 a[LDA*N] = {
	   { 9.14,  0.00}, { 0.00,  0.00}, { 0.00,  0.00}, { 0.00,  0.00},
	   {-4.37,  9.22}, {-3.35,  0.00}, { 0.00,  0.00}, { 0.00,  0.00},
	   {-1.98,  1.72}, { 2.25,  9.51}, {-4.82,  0.00}, { 0.00,  0.00},
	   {-8.96,  9.50}, { 2.57, -2.40}, {-3.24, -2.04}, { 8.44,  0.00}
	};
	/* Executable statements */
	printf( "LAPACKE_zheev (row-major, high-level) Example Program Results\n" );
	/* Solve eigenproblem */
	info = LAPACKE_zheev( LAPACK_ROW_MAJOR, 'V', 'L', n, a, lda, w );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	/* Print eigenvalues */
	print_rmatrix( "Eigenvalues", 1, n, w, 1 );
	/* Print eigenvectors */
	print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
	exit( 0 );
} /* End of LAPACKE_zheev Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%6.2f,%6.2f)", a[i*lda+j].real, a[i*lda+j].imag );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
		printf( "\n" );
	}
}
