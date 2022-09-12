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
   SSYEVR Example.
   ==============

   Program computes the smallest eigenvalues and the corresponding
   eigenvectors of a real symmetric matrix A using the Relatively Robust
   Representations, where A is:

     0.67  -0.20   0.19  -1.06   0.46
    -0.20   3.82  -0.13   1.06  -0.48
     0.19  -0.13   3.27   0.11   1.10
    -1.06   1.06   0.11   5.86  -0.98
     0.46  -0.48   1.10  -0.98   3.54

   Description.
   ============

   The routine computes selected eigenvalues and, optionally, eigenvectors of
   an n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.
   Eigenvalues and eigenvectors can be selected by specifying either a range
   of values or a range of indices for the desired eigenvalues.

   Example Program Results.
   ========================

 SSYEVR Example Program Results

 The total number of eigenvalues found: 3

 Selected eigenvalues
   0.43   2.14   3.37

 Selected eigenvectors (stored columnwise)
  -0.98  -0.01  -0.08
   0.01   0.02  -0.93
   0.04  -0.69  -0.07
  -0.18   0.19   0.31
   0.07   0.69  -0.13
*/
#include <stdlib.h>
#include <stdio.h>
#include <mkl.h>

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, int m, int n, float* a, int lda );

/* Parameters */
#define N 5
#define NSELECT 3
#define LDA N
#define LDZ N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, il, iu, m, lda = LDA, ldz = LDZ, info, lwork, liwork;
	float abstol, vl, vu;
	MKL_INT iwkopt;
	MKL_INT* iwork;
	float wkopt;
	float* work;
	/* Local arrays */
	MKL_INT isuppz[N];
	float w[N], z[LDZ*NSELECT];
	float a[LDA*N] = {
	    0.67f,  0.00f,  0.00f,  0.00f,  0.00f,
	   -0.20f,  3.82f,  0.00f,  0.00f,  0.00f,
	    0.19f, -0.13f,  3.27f,  0.00f,  0.00f,
	   -1.06f,  1.06f,  0.11f,  5.86f,  0.00f,
	    0.46f, -0.48f,  1.10f, -0.98f,  3.54f
	};
	/* Executable statements */
	printf( " SSYEVR Example Program Results\n" );
	/* Negative abstol means using the default value */
	abstol = -1.0;
	/* Set il, iu to compute NSELECT smallest eigenvalues */
	il = 1;
	iu = NSELECT;
	/* Query and allocate the optimal workspace */
	lwork = -1;
	liwork = -1;
	ssyevr( "Vectors", "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &m, w, z, &ldz, isuppz, &wkopt, &lwork, &iwkopt, &liwork,
			&info );
	lwork = (MKL_INT)wkopt;
	work = (float*)malloc( lwork*sizeof(float) );
	liwork = iwkopt;
	iwork = (MKL_INT*)malloc( liwork*sizeof(MKL_INT) );
	/* Solve eigenproblem */
	ssyevr( "Vectors", "Indices", "Upper", &n, a, &lda, &vl, &vu, &il, &iu,
			&abstol, &m, w, z, &ldz, isuppz, work, &lwork, iwork, &liwork,
			&info );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	/* Print the number of eigenvalues found */
	printf( "\n The total number of eigenvalues found:%2i\n", m );
	/* Print eigenvalues */
	print_matrix( "Selected eigenvalues", 1, m, w, 1 );
	/* Print eigenvectors */
	print_matrix( "Selected eigenvectors (stored columnwise)", n, m, z, ldz );
	/* Free workspace */
	free( (void*)iwork );
	free( (void*)work );
	exit( 0 );
} /* End of SSYEVR Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, int m, int n, float* a, int lda ) {
	int i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}
