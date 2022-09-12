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
   SPOSV Example.
   ==============
 
   The program computes the solution to the system of linear
   equations with a symmetric positive-definite matrix A and multiple
   right-hand sides B, where A is the coefficient matrix:
 
     3.14   0.17  -0.90   1.65  -0.72
     0.17   0.79   0.83  -0.65   0.28
    -0.90   0.83   4.53  -3.70   1.60
     1.65  -0.65  -3.70   5.32  -1.37
    -0.72   0.28   1.60  -1.37   1.98

   and B is the right-hand side matrix:
 
    -7.29   6.11   0.59
     9.25   2.90   8.88
     5.99  -5.05   7.57
    -1.94  -3.80   5.57
    -8.30   9.66  -1.67
 
   Description.
   ============
 
   The routine solves for X the real system of linear equations 
   A*X = B, where A is an n-by-n symmetric positive-definite 
   matrix, the columns of matrix B are individual right-hand sides, 
   and the columns of X are the corresponding solutions.

   The Cholesky decomposition is used to factor A as
   A = UT*U, if uplo = 'U' or A = L*LT, if uplo = 'L',
   where U is an upper triangular matrix and L is a lower triangular matrix.
   The factored form of A is then used to solve the system of equations A*X = B.

   Example Program Results.
   ========================
 
 SPOSV Example Program Results

 Solution
  -6.02   3.95  -3.14
  15.62   4.32  13.05
   3.02  -8.25   4.91
   3.25  -4.83   6.11
  -8.78   9.04  -3.57

 Details of Cholesky factorization
   1.77   0.10  -0.51   0.93  -0.41
   0.00   0.88   0.99  -0.84   0.36
   0.00   0.00   1.81  -1.32   0.57
   0.00   0.00   0.00   1.42   0.05
   0.00   0.00   0.00   0.00   1.16
*/
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, float* a, int lda );

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
	float a[LDA*N] = {
	    3.14f,  0.00f,  0.00f,  0.00f,  0.00f,
	    0.17f,  0.79f,  0.00f,  0.00f,  0.00f,
	   -0.90f,  0.83f,  4.53f,  0.00f,  0.00f,
	    1.65f, -0.65f, -3.70f,  5.32f,  0.00f,
	   -0.72f,  0.28f,  1.60f, -1.37f,  1.98f
	};
	float b[LDB*NRHS] = {
	   -7.29f,  9.25f,  5.99f, -1.94f, -8.30f,
	    6.11f,  2.90f, -5.05f, -3.80f,  9.66f,
	    0.59f,  8.88f,  7.57f,  5.57f, -1.67f
	};
	/* Executable statements */
	printf( " SPOSV Example Program Results\n" );
	/* Solve the equations A*X = B */
	sposv( "Upper", &n, &nrhs, a, &lda, b, &ldb, &info );
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
} /* End of SPOSV Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, float* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}
