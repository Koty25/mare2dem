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
   LAPACKE_zheevx Example.
   =======================

   Program computes eigenvalues specified by a selected range of values
   and corresponding eigenvectors of a complex Hermitian matrix A:

   (  6.51,  0.00) ( -5.92,  9.53) ( -2.46,  2.91) (  8.84,  3.21)
   ( -5.92, -9.53) ( -1.73,  0.00) (  6.50,  2.09) (  1.32,  8.81)
   ( -2.46, -2.91) (  6.50, -2.09) (  6.90,  0.00) ( -0.59,  2.47)
   (  8.84, -3.21) (  1.32, -8.81) ( -0.59, -2.47) ( -2.85,  0.00)

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

 LAPACKE_zheevx (column-major, high-level) Example Program Results

 The total number of eigenvalues found: 3

 Selected eigenvalues
   0.09   9.53  18.75

 Selected eigenvectors (stored columnwise)
 (  0.18,  0.00) ( -0.54,  0.00) (  0.67,  0.00)
 ( -0.40, -0.31) ( -0.21, -0.17) ( -0.30, -0.43)
 (  0.60,  0.40) ( -0.35, -0.28) ( -0.39, -0.34)
 ( -0.34,  0.26) ( -0.57,  0.35) (  0.05,  0.05)
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
#define LDZ N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, lda = LDA, ldz = LDZ, il, iu, m, info;
	double abstol, vl, vu;
	/* Local arrays */
	MKL_INT ifail[N];
	double w[N];
	MKL_Complex16 z[LDZ*N];
	MKL_Complex16 a[LDA*N] = {
	   { 6.51,  0.00}, {-5.92, -9.53}, {-2.46, -2.91}, { 8.84, -3.21},
	   { 0.00,  0.00}, {-1.73,  0.00}, { 6.50, -2.09}, { 1.32, -8.81},
	   { 0.00,  0.00}, { 0.00,  0.00}, { 6.90,  0.00}, {-0.59, -2.47},
	   { 0.00,  0.00}, { 0.00,  0.00}, { 0.00,  0.00}, {-2.85,  0.00}
	};
	/* Executable statements */
	printf( "LAPACKE_zheevx (column-major, high-level) Example Program Results\n" );
	/* Negative abstol means using the default value */
	abstol = -1.0;
	/* Set VL, VU to compute eigenvalues in half-open (VL,VU] interval */
	vl = 0.0;
	vu = 100.0;
	/* Solve eigenproblem */
	info = LAPACKE_zheevx( LAPACK_COL_MAJOR, 'V', 'V', 'L', n, a, lda,
			vl, vu, il, iu, abstol, &m, w, z, ldz, ifail );
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
	exit( 0 );
} /* End of LAPACKE_zheevx Example */

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
