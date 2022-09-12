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
   LAPACKE_cgels Example.
   ======================

   Program computes the minimum norm solution to the underdetermined linear
   system A*X = B with full rank matrix A using LQ factorization,
   where A is the coefficient matrix:

   ( -4.20, -3.44) ( -3.35,  1.52) (  1.73,  8.85) (  2.35,  0.34)
   ( -5.43, -8.81) ( -4.53, -8.47) (  5.93,  3.75) ( -3.75, -5.66)
   ( -5.56,  3.39) (  2.90, -9.22) (  8.03,  9.37) (  5.69, -0.47)

   and B is the right-hand side matrix:

   ( -7.02,  4.80) (  3.88, -2.59)
   (  0.62, -2.40) (  1.57,  3.24)
   (  3.10, -2.19) ( -6.93, -5.99)

   Description.
   ============

   The routine solves overdetermined or underdetermined complex linear systems
   involving an m-by-n matrix A, or its transpose, using a QR or LQ
   factorization of A. It is assumed that A has full rank.

   Several right hand side vectors b and solution vectors x can be handled
   in a single call; they are stored as the columns of the m-by-nrhs right
   hand side matrix B and the n-by-nrhs solution matrix X.

   Example Program Results.
   ========================

 LAPACKE_cgels (column-major, high-level) Example Program Results

 Minimum norm solution
 ( -0.25, -0.04) ( -0.21,  0.42)
 (  0.99,  0.27) ( -0.21, -0.43)
 (  0.25,  0.43) ( -0.24, -0.13)
 ( -0.32,  0.14) ( -0.23, -0.09)

 Details of LQ factorization
 ( 11.40,  0.00) (  0.18, -0.14) ( -0.23, -0.52) ( -0.15,  0.01)
 (  7.73, -0.39) ( 15.32,  0.00) ( -0.22,  0.42) (  0.45,  0.17)
 (  8.60, -5.68) (  3.96,  6.46) ( 12.54,  0.00) ( -0.02, -0.47)
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda );

/* Parameters */
#define M 3
#define N 4
#define NRHS 2
#define LDA M
#define LDB N

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
	/* Local arrays */
	MKL_Complex8 a[LDA*N] = {
	   {-4.20f, -3.44f}, {-5.43f, -8.81f}, {-5.56f,  3.39f},
	   {-3.35f,  1.52f}, {-4.53f, -8.47f}, { 2.90f, -9.22f},
	   { 1.73f,  8.85f}, { 5.93f,  3.75f}, { 8.03f,  9.37f},
	   { 2.35f,  0.34f}, {-3.75f, -5.66f}, { 5.69f, -0.47f}
	};
	MKL_Complex8 b[LDB*NRHS] = {
	   {-7.02f,  4.80f}, { 0.62f, -2.40f}, { 3.10f, -2.19f}, { 0.00f,  0.00f},
	   { 3.88f, -2.59f}, { 1.57f,  3.24f}, {-6.93f, -5.99f}, { 0.00f,  0.00f}
	};
	/* Executable statements */
	printf( "LAPACKE_cgels (column-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */
	info = LAPACKE_cgels( LAPACK_COL_MAJOR, 'N', m, n, nrhs, a, lda, b, ldb );
	/* Check for the full rank */
	if( info > 0 ) {
		printf( "The diagonal element %i of the triangular factor ", info );
		printf( "of A is zero, so that A does not have full rank;\n" );
		printf( "the minimum norm solution could not be computed.\n" );
		exit( 1 );
	}
	/* Print minimum norm solution */
	print_matrix( "Minimum norm solution", n, nrhs, b, ldb );
	/* Print details of LQ factorization */
	print_matrix( "Details of LQ factorization", m, n, a, lda );
	exit( 0 );
} /* End of LAPACKE_cgels Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%6.2f,%6.2f)", a[i+j*lda].real, a[i+j*lda].imag );
		printf( "\n" );
	}
}
