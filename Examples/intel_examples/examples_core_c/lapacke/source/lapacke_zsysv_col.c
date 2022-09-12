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
   LAPACKE_zsysv Example.
   ======================
 
   The program computes the solution to the system of linear equations
   with a complex symmetric matrix A and multiple right-hand sides B,
   where A is the coefficient matrix:

   (  9.99, -4.73) ( -5.68, -0.80) ( -8.94,  1.32) ( -9.42,  2.05)
   ( -5.68, -0.80) ( -8.01,  4.61) (  1.64, -6.29) (  6.79, -2.17)
   ( -8.94,  1.32) (  1.64, -6.29) (  9.04,  3.96) ( -4.51, -7.54)
   ( -9.42,  2.05) (  6.79, -2.17) ( -4.51, -7.54) (  0.40,  4.06)

   and B is the right-hand side matrix:
 
   (  5.71, -1.20) (  2.84, -0.18)
   ( -7.70,  6.47) ( -8.29, -1.72)
   (  3.77, -7.40) ( -4.28, -8.25)
   ( -3.78,  0.33) ( -2.70, -0.39)
 
   Description.
   ============
 
   The routine solves for X the complex system of linear equations A*X = B,
   where A is an n-by-n symmetric matrix, the columns of matrix B are
   individual right-hand sides, and the columns of X are the corresponding
   solutions.

   The diagonal pivoting method is used to factor A as A = U*D*UT or
   A = L*D*LT , where U (or L) is a product of permutation and unit upper
   (lower) triangular matrices, and D is symmetric and block diagonal with
   1-by-1 and 2-by-2 diagonal blocks.

   The factored form of A is then used to solve the system of equations A*X = B.

   Example Program Results.
   ========================
 
 LAPACKE_zsysv (column-major, high-level) Example Program Results

 Solution
 (  0.13,  0.13) (  0.63,  0.34)
 (  0.32, -0.07) (  0.61,  0.21)
 ( -0.26, -0.44) ( -0.01, -0.10)
 ( -0.40,  0.51) (  0.21,  0.02)

 Details of factorization
 (-16.42,  1.69) ( -0.53,  0.35) (  0.36,  0.41) ( -0.78,  0.49)
 (  0.00,  0.00) (  3.69,  0.64) (-16.58, -1.61) ( -0.10, -0.65)
 (  0.00,  0.00) (  0.00,  0.00) (  1.02, -3.74) ( -0.73, -0.52)
 (  0.00,  0.00) (  0.00,  0.00) (  0.00,  0.00) (  9.04,  3.96)

 Pivot indices
      1     -1     -1      3
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda );
extern void print_int_vector( char* desc, MKL_INT n, MKL_INT* a );

/* Parameters */
#define N 4
#define NRHS 2
#define LDA N
#define LDB N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
	/* Local arrays */
	MKL_INT ipiv[N];
	MKL_Complex16 a[LDA*N] = {
	   { 9.99, -4.73}, { 0.00,  0.00}, { 0.00,  0.00}, { 0.00,  0.00},
	   {-5.68, -0.80}, {-8.01,  4.61}, { 0.00,  0.00}, { 0.00,  0.00},
	   {-8.94,  1.32}, { 1.64, -6.29}, { 9.04,  3.96}, { 0.00,  0.00},
	   {-9.42,  2.05}, { 6.79, -2.17}, {-4.51, -7.54}, { 0.40,  4.06}
	};
	MKL_Complex16 b[LDB*NRHS] = {
	   { 5.71, -1.20}, {-7.70,  6.47}, { 3.77, -7.40}, {-3.78,  0.33},
	   { 2.84, -0.18}, {-8.29, -1.72}, {-4.28, -8.25}, {-2.70, -0.39}
	};
	/* Executable statements */
	printf( "LAPACKE_zsysv (column-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */
	info = LAPACKE_zsysv( LAPACK_COL_MAJOR, 'U', n, nrhs, a, lda, ipiv,
			b, ldb );
	/* Check for the exact singularity */
	if( info > 0 ) {
		printf( "The element of the diagonal factor " );
		printf( "D(%i,%i) is zero, so that D is singular;\n", info, info );
		printf( "the solution could not be computed.\n" );
		exit( 1 );
	}
	/* Print solution */
	print_matrix( "Solution", n, nrhs, b, ldb );
	/* Print details of factorization */
	print_matrix( "Details of factorization", n, n, a, lda );
	/* Print pivot indices */
	print_int_vector( "Pivot indices", n, ipiv );
	exit( 0 );
} /* End of LAPACKE_zsysv Example */

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

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, MKL_INT n, MKL_INT* a ) {
	MKL_INT j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
}
