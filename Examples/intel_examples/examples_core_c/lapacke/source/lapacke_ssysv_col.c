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
   LAPACKE_ssysv Example.
   ======================
 
   The program computes the solution to the system of linear equations
   with a real symmetric matrix A and multiple right-hand sides B,
   where A is the coefficient matrix:

    -5.86   3.99  -5.93  -2.82   7.69
     3.99   4.46   2.58   4.42   4.61
    -5.93   2.58  -8.52   8.57   7.69
    -2.82   4.42   8.57   3.72   8.07
     7.69   4.61   7.69   8.07   9.83

   and B is the right-hand side matrix:
 
     1.32  -6.33  -8.77
     2.22   1.69  -8.33
     0.12  -1.56   9.54
    -6.41  -9.49   9.56
     6.33  -3.67   7.48
 
   Description.
   ============
 
   The routine solves for X the real system of linear equations A*X = B,
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
 
 LAPACKE_ssysv (column-major, high-level) Example Program Results

 Solution
   1.17   0.52  -0.86
  -0.71   1.05  -4.90
  -0.63  -0.52   0.99
  -0.33   0.43   1.22
   0.83  -1.22   1.96

 Details of factorization
  -5.86   0.00   0.00   0.00   0.00
  -0.68   7.18   0.00   0.00   0.00
   1.01  -0.20  -2.82   0.00   0.00
   0.48   0.35  11.93   4.21   0.00
  -1.31   1.37   0.02   0.16   6.22

 Pivot indices
      1      2     -4     -4      5
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda );
extern void print_int_vector( char* desc, MKL_INT n, MKL_INT* a );

/* Parameters */
#define N 5
#define NRHS 3
#define LDA N
#define LDB N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, nrhs = NRHS, lda = LDA, ldb = LDB, info;
	/* Local arrays */
	MKL_INT ipiv[N];
	float a[LDA*N] = {
	   -5.86f,  3.99f, -5.93f, -2.82f,  7.69f,
	    0.00f,  4.46f,  2.58f,  4.42f,  4.61f,
	    0.00f,  0.00f, -8.52f,  8.57f,  7.69f,
	    0.00f,  0.00f,  0.00f,  3.72f,  8.07f,
	    0.00f,  0.00f,  0.00f,  0.00f,  9.83f
	};
	float b[LDB*NRHS] = {
	    1.32f,  2.22f,  0.12f, -6.41f,  6.33f,
	   -6.33f,  1.69f, -1.56f, -9.49f, -3.67f,
	   -8.77f, -8.33f,  9.54f,  9.56f,  7.48f
	};
	/* Executable statements */
	printf( "LAPACKE_ssysv (column-major, high-level) Example Program Results\n" );
	/* Solve the equations A*X = B */
	info = LAPACKE_ssysv( LAPACK_COL_MAJOR, 'L', n, nrhs, a, lda, ipiv,
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
} /* End of LAPACKE_ssysv Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
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
