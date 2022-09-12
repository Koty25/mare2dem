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
   LAPACKE_zgesvd Example.
   =======================

   Program computes the singular value decomposition of a general
   rectangular complex matrix A:

   (  5.91, -5.69) (  7.09,  2.72) (  7.78, -4.06) ( -0.79, -7.21)
   ( -3.15, -4.08) ( -1.89,  3.27) (  4.57, -2.07) ( -3.88, -3.30)
   ( -4.89,  4.20) (  4.10, -6.70) (  3.28, -3.84) (  3.84,  1.19)

   Description.
   ============

   The routine computes the singular value decomposition (SVD) of a complex
   m-by-n matrix A, optionally computing the left and/or right singular
   vectors. The SVD is written as

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

 LAPACKE_zgesvd (column-major, high-level) Example Program Results

 Singular values
  17.63  11.61   6.78

 Left singular vectors (stored columnwise)
 ( -0.86,  0.00) (  0.40,  0.00) (  0.32,  0.00)
 ( -0.35,  0.13) ( -0.24, -0.21) ( -0.63,  0.60)
 (  0.15,  0.32) (  0.61,  0.61) ( -0.36,  0.10)

 Right singular vectors (stored rowwise)
 ( -0.22,  0.51) ( -0.37, -0.32) ( -0.53,  0.11) (  0.15,  0.38)
 (  0.31,  0.31) (  0.09, -0.57) (  0.18, -0.39) (  0.38, -0.39)
 (  0.53,  0.24) (  0.49,  0.28) ( -0.47, -0.25) ( -0.15,  0.19)
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

#define min(a,b) ((a)>(b)?(b):(a))

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda );
extern void print_rmatrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
#define M 3
#define N 4
#define LDA M
#define LDU M
#define LDVT N

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info;
	/* Local arrays */
	double s[M];
	double superb[min(M,N)-1];
	MKL_Complex16 u[LDU*M], vt[LDVT*N];
	MKL_Complex16 a[LDA*N] = {
	   { 5.91, -5.69}, {-3.15, -4.08}, {-4.89,  4.20},
	   { 7.09,  2.72}, {-1.89,  3.27}, { 4.10, -6.70},
	   { 7.78, -4.06}, { 4.57, -2.07}, { 3.28, -3.84},
	   {-0.79, -7.21}, {-3.88, -3.30}, { 3.84,  1.19}
	};
	/* Executable statements */
	printf( "LAPACKE_zgesvd (column-major, high-level) Example Program Results\n" );
	/* Compute SVD */
	info = LAPACKE_zgesvd( LAPACK_COL_MAJOR, 'A', 'A', m, n, a, lda, s,
         u, ldu, vt, ldvt, superb );
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
} /* End of LAPACKE_zgesvd Example */

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
