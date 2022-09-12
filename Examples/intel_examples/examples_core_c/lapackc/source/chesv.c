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
   CHESV Example.
   ==============
 
   The program computes the solution to the system of linear equations
   with a Hermitian matrix A and multiple right-hand sides B,
   where A is the coefficient matrix:

   ( -2.90,  0.00) (  0.31,  4.46) (  9.66, -5.66) ( -2.28,  2.14)
   (  0.31, -4.46) ( -7.93,  0.00) (  9.55, -4.62) ( -3.51,  3.11)
   (  9.66,  5.66) (  9.55,  4.62) (  0.30,  0.00) (  9.33, -9.66)
   ( -2.28, -2.14) ( -3.51, -3.11) (  9.33,  9.66) (  2.40,  0.00)

   and B is the right-hand side matrix:
 
   ( -5.69, -8.21) ( -2.83,  6.46)
   ( -3.57,  1.99) ( -7.64,  1.10)
   (  8.42, -9.83) ( -2.33, -4.23)
   ( -5.00,  3.85) (  6.48, -3.81)
 
   Description.
   ============
 
   The routine solves for X the complex system of linear equations A*X = B,
   where A is an n-by-n Hermitian matrix, the columns of matrix B are
   individual right-hand sides, and the columns of X are the corresponding
   solutions.

   The diagonal pivoting method is used to factor A as A = U*D*UH or
   A = L*D*LH, where U (or L) is a product of permutation and unit upper
   (lower) triangular matrices, and D is Hermitian and block diagonal with
   1-by-1 and 2-by-2 diagonal blocks.

   The factored form of A is then used to solve the system of equations A*X = B.

   Example Program Results.
   ========================
 
 CHESV Example Program Results

 Solution
 (  0.22, -0.95) ( -1.13,  0.18)
 ( -1.42, -1.30) (  0.70,  1.13)
 ( -0.65, -0.40) (  0.04,  0.07)
 ( -0.48,  1.35) (  1.15, -0.27)

 Details of factorization
 (  3.17,  0.00) (  7.32,  3.28) ( -0.36,  0.06) (  0.20, -0.82)
 (  0.00,  0.00) (  0.03,  0.00) ( -0.48,  0.03) (  0.25, -0.76)
 (  0.00,  0.00) (  0.00,  0.00) (  0.30,  0.00) (  9.33, -9.66)
 (  0.00,  0.00) (  0.00,  0.00) (  0.00,  0.00) (  2.40,  0.00)

 Pivot indices
     -1     -1     -3     -3
*/
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>

/* Complex datatype */
struct _fcomplex { float re, im; };
typedef struct _fcomplex fcomplex;

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, fcomplex* a, int lda );
extern void print_int_vector( char* desc, int n, int* a );

/* Parameters */
#define N 4
#define NRHS 2
#define LDA N
#define LDB N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info, lwork;
	fcomplex wkopt;
	fcomplex* work;
	/* Local arrays */
	MKL_INT ipiv[N];
	fcomplex a[LDA*N] = {
	   {-2.90f,  0.00f}, { 0.00f,  0.00f}, { 0.00f,  0.00f}, { 0.00f,  0.00f},
	   { 0.31f,  4.46f}, {-7.93f,  0.00f}, { 0.00f,  0.00f}, { 0.00f,  0.00f},
	   { 9.66f, -5.66f}, { 9.55f, -4.62f}, { 0.30f,  0.00f}, { 0.00f,  0.00f},
	   {-2.28f,  2.14f}, {-3.51f,  3.11f}, { 9.33f, -9.66f}, { 2.40f,  0.00f}
	};
	fcomplex b[LDB*NRHS] = {
	   {-5.69f, -8.21f}, {-3.57f,  1.99f}, { 8.42f, -9.83f}, {-5.00f,  3.85f},
	   {-2.83f,  6.46f}, {-7.64f,  1.10f}, {-2.33f, -4.23f}, { 6.48f, -3.81f}
	};
	/* Executable statements */
	printf( " CHESV Example Program Results\n" );
	/* Query and allocate the optimal workspace */
	lwork = -1;
	chesv( "Upper", &n, &nrhs, a, &lda, ipiv, b, &ldb, &wkopt, &lwork, &info );
	lwork = (MKL_INT)wkopt.re;
	work = (fcomplex*)malloc( lwork*sizeof(fcomplex) );
	/* Solve the equations A*X = B */
	chesv( "Upper", &n, &nrhs, a, &lda, ipiv, b, &ldb, work, &lwork, &info );
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
	/* Free workspace */
	free( (void*)work );
	exit( 0 );
} /* End of CHESV Example */

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

/* Auxiliary routine: printing a vector of integers */
void print_int_vector( char* desc, int n, int* a ) {
	int j;
	printf( "\n %s\n", desc );
	for( j = 0; j < n; j++ ) printf( " %6i", a[j] );
	printf( "\n" );
}
