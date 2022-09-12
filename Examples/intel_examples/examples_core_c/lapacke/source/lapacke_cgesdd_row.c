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
   LAPACKE_cgesdd Example.
   =======================

   Program computes the singular value decomposition of a general
   rectangular complex matrix A using a divide and conquer method, where A is:

   ( -5.40,  7.40) (  6.00,  6.38) (  9.91,  0.16) ( -5.28, -4.16)
   (  1.09,  1.55) (  2.60,  0.07) (  3.98, -5.26) (  2.03,  1.11)
   (  9.88,  1.91) (  4.92,  6.31) ( -2.11,  7.39) ( -9.81, -8.98)

   Description.
   ============

   The routine computes the singular value decomposition (SVD) of a complex
   m-by-n matrix A, optionally computing the left and/or right singular
   vectors. If singular vectors are desired, it uses a divide and conquer
   algorithm. The SVD is written as

   A = U*SIGMA*VH

   where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
   diagonal elements, U is an m-by-m unitary matrix and VH (V conjugate
   transposed) is an n-by-n unitary matrix. The diagonal elements of SIGMA
   are the singular values of A; they are real and non-negative, and are
   returned in descending order. The first min(m, n) columns of U and V are
   the left and right singular vectors of A.

   Note that the routine returns VH, not V.

   Example Program Results.
   ========================

 LAPACKE_cgesdd (row-major, high-level) Example Program Results

 Singular values
  21.76  16.60   3.97

 Left singular vectors (stored columnwise)
 (  0.55,  0.00) (  0.76,  0.00) ( -0.34,  0.00)
 ( -0.04, -0.15) (  0.27, -0.23) (  0.55, -0.74)
 (  0.81,  0.12) ( -0.52, -0.14) (  0.13, -0.11)

 Right singular vectors (stored rowwise)
 (  0.23,  0.21) (  0.37,  0.39) (  0.24,  0.33) ( -0.56, -0.37)
 ( -0.58,  0.40) (  0.11,  0.17) (  0.60, -0.27) (  0.16,  0.06)
 (  0.60,  0.12) ( -0.19,  0.30) (  0.39,  0.20) (  0.45,  0.31)
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda );
extern void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda );

/* Parameters */
#define M 3
#define N 4
#define LDA N
#define LDU M
#define LDVT N

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info;
	/* Local arrays */
	float s[M];
	MKL_Complex8 u[LDU*M], vt[LDVT*N];
	MKL_Complex8 a[LDA*M] = {
	   {-5.40f,  7.40f}, { 6.00f,  6.38f}, { 9.91f,  0.16f}, {-5.28f, -4.16f},
	   { 1.09f,  1.55f}, { 2.60f,  0.07f}, { 3.98f, -5.26f}, { 2.03f,  1.11f},
	   { 9.88f,  1.91f}, { 4.92f,  6.31f}, {-2.11f,  7.39f}, {-9.81f, -8.98f}
	};
	/* Executable statements */
	printf( "LAPACKE_cgesdd (row-major, high-level) Example Program Results\n" );
	/* Compute SVD */
	info = LAPACKE_cgesdd( LAPACK_ROW_MAJOR, 'S', m, n, a, lda, s, u, 
         ldu, vt, ldvt );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
	}
	/* Print singular values */
	print_rmatrix( "Singular values", 1, m, s, 1 );
	/* Print left singular vectors */
	print_matrix( "Left singular vectors (stored columnwise)", m, m, u, ldu );
	/* Print right singular vectors */
	print_matrix( "Right singular vectors (stored rowwise)", m, n, vt, ldvt );
	exit( 0 );
} /* End of LAPACKE_cgesdd Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ )
			printf( " (%6.2f,%6.2f)", a[i*lda+j].real, a[i*lda+j].imag );
		printf( "\n" );
	}
}

/* Auxiliary routine: printing a real matrix */
void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
		printf( "\n" );
	}
}
