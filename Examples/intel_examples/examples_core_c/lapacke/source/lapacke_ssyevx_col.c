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
   LAPACKE_ssyevx Example.
   =======================

   Program computes the smallest eigenvalues and the corresponding
   eigenvectors of a real symmetric matrix A:

     6.29  -0.39   0.61   1.18  -0.08
    -0.39   7.19   0.81   1.19  -0.08
     0.61   0.81   5.48  -3.13   0.22
     1.18   1.19  -3.13   3.79  -0.26
    -0.08  -0.08   0.22  -0.26   0.83

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

 LAPACKE_ssyevx (column-major, high-level) Example Program Results

 The total number of eigenvalues found: 3

 Selected eigenvalues
   0.71   0.82   6.58

 Selected eigenvectors (stored columnwise)
   0.22   0.09  -0.95
   0.21   0.08  -0.04
  -0.52  -0.22  -0.29
  -0.73  -0.21  -0.09
  -0.32   0.94   0.01
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda );

/* Parameters */
#define N 5
#define NSELECT 3
#define LDA N
#define LDZ N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, il, iu, m, lda = LDA, ldz = LDZ, info;
	float abstol, vl = 0.0f, vu = 0.0f;
	/* Local arrays */
	MKL_INT ifail[N];
	float w[N], z[LDZ*NSELECT];
	float a[LDA*N] = {
	    6.29f,  0.00f,  0.00f,  0.00f,  0.00f,
	   -0.39f,  7.19f,  0.00f,  0.00f,  0.00f,
	    0.61f,  0.81f,  5.48f,  0.00f,  0.00f,
	    1.18f,  1.19f, -3.13f,  3.79f,  0.00f,
	   -0.08f, -0.08f,  0.22f, -0.26f,  0.83f
	};
	/* Executable statements */
	printf( "LAPACKE_ssyevx (column-major, high-level) Example Program Results\n" );
	/* Negative abstol means using the default value */
	abstol = -1.0;
	/* Set il, iu to compute NSELECT smallest eigenvalues */
	il = 1;
	iu = NSELECT;
	/* Solve eigenproblem */
	info = LAPACKE_ssyevx( LAPACK_COL_MAJOR, 'V', 'I', 'U', n, a, lda,
			vl, vu, il, iu, abstol, &m, w, z, ldz, ifail );
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
	exit( 0 );
} /* End of LAPACKE_ssyevx Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}
