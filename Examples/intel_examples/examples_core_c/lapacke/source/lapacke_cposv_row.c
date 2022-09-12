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
   LAPACKE_cposv Example.
   ======================
 
   The program computes the solution to the system of linear
   equations with a Hermitian positive-definite matrix A and multiple
   right-hand sides B, where A is the coefficient matrix:
 
   (  5.96,  0.00) (  0.40, -1.19) ( -0.83, -0.48) ( -0.57,  0.40)
   (  0.40,  1.19) (  7.95,  0.00) (  0.33,  0.09) (  0.22,  0.74)
   ( -0.83,  0.48) (  0.33, -0.09) (  4.43,  0.00) ( -1.09,  0.32)
   ( -0.57, -0.40) (  0.22, -0.74) ( -1.09, -0.32) (  3.46,  0.00)

   and B is the right-hand side matrix:
 
   ( -2.94,  5.79) (  8.44,  3.07)
   (  8.12, -9.12) (  1.00, -4.62)
   (  9.09, -5.03) (  3.64, -2.33)
   (  7.36,  6.77) (  8.04,  2.87)
 
   Description.
   ============
 
   The routine solves for X the complex system of linear equations 
   A*X = B, where A is an n-by-n Hermitian positive-definite 
   matrix, the columns of matrix B are individual right-hand sides, 
   and the columns of X are the corresponding solutions.

   The Cholesky decomposition is used to factor A as
   A = UH*U, if uplo = 'U' or A = L*LH, if uplo = 'L',
   where U is an upper triangular matrix and L is a lower triangular matrix.
   The factored form of A is then used to solve the system of equations A*X = B.

   Example Program Results.
   ========================
 
 LAPACKE_cposv (row-major, high-level) Example Program Results

 Solution
 (  0.80,  1.62) (  2.52,  0.61)
 (  1.26, -1.78) (  0.01, -1.38)
 (  3.38, -0.29) (  2.42, -0.52)
 (  3.46,  2.92) (  3.77,  1.37)

 Details of Cholesky factorization
 (  2.44,  0.00) (  0.00,  0.00) (  0.00,  0.00) (  0.00,  0.00)
 (  0.16,  0.49) (  2.77,  0.00) (  0.00,  0.00) (  0.00,  0.00)
 ( -0.34,  0.20) (  0.10, -0.10) (  2.06,  0.00) (  0.00,  0.00)
 ( -0.23, -0.16) (  0.12, -0.30) ( -0.57, -0.20) (  1.71,  0.00)
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda );

/* Parameters */
#define N 4
#define NRHS 2
#define LDA N
#define LDB NRHS

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
	/* Local arrays */
	MKL_Complex8 a[LDA*N] = {
	   { 5.96f,  0.00f}, { 0.00f,  0.00f}, { 0.00f,  0.00f}, { 0.00f,  0.00f},
	   { 0.40f,  1.19f}, { 7.95f,  0.00f}, { 0.00f,  0.00f}, { 0.00f,  0.00f},
	   {-0.83f,  0.48f}, { 0.33f, -0.09f}, { 4.43f,  0.00f}, { 0.00f,  0.00f},
	   {-0.57f, -0.40f}, { 0.22f, -0.74f}, {-1.09f, -0.32f}, { 3.46f,  0.00f}
	};
	MKL_Complex8 b[LDB*N] = {
	   {-2.94f,  5.79f}, { 8.44f,  3.07f},
	   { 8.12f, -9.12f}, { 1.00f, -4.62f},
	   { 9.09f, -5.03f}, { 3.64f, -2.33f},
	   { 7.36f,  6.77f}, { 8.04f,  2.87f}
	};
	/* Executable statements */
	printf( "LAPACKE_cposv (row-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */
	info = LAPACKE_cposv( LAPACK_ROW_MAJOR, 'L', n, nrhs, a, lda, b, ldb );
	/* Check for the positive definiteness */
	if( info > 0 ) {
		printf( "The leading minor of order %i is not positive ", info );
		printf( "definite;\nthe solution could not be computed.\n" );
		exit( 1 );
	}
	/* Print solution */
	print_matrix( "Solution", n, nrhs, b, ldb );
	/* Print details of Cholesky factorization */
	print_matrix( "Details of Cholesky factorization", n, n, a, lda );
	exit( 0 );
} /* End of LAPACKE_cposv Example */

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
