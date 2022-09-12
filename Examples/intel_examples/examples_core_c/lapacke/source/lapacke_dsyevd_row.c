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
   LAPACKE_dsyevd Example.
   =======================

   Program computes all eigenvalues and eigenvectors of a real symmetric
   matrix A using divide and conquer algorithm, where A is:

     6.39   0.13  -8.23   5.71  -3.18
     0.13   8.37  -4.46  -6.10   7.21
    -8.23  -4.46  -9.58  -9.25  -7.42
     5.71  -6.10  -9.25   3.72   8.54
    -3.18   7.21  -7.42   8.54   2.51

   Description.
   ============

   The routine computes all eigenvalues and, optionally, eigenvectors of an
   n-by-n real symmetric matrix A. The eigenvector v(j) of A satisfies

   A*v(j) = lambda(j)*v(j)

   where lambda(j) is its eigenvalue. The computed eigenvectors are
   orthonormal.
   If the eigenvectors are requested, then this routine uses a divide and
   conquer algorithm to compute eigenvalues and eigenvectors.

   Example Program Results.
   ========================

 LAPACKE_dsyevd (row-major, high-level) Example Program Results

 Eigenvalues
 -17.44 -11.96   6.72  14.25  19.84

 Eigenvectors (stored columnwise)
  -0.26   0.31  -0.74   0.33   0.42
  -0.17  -0.39  -0.38  -0.80   0.16
  -0.89   0.04   0.09   0.03  -0.45
  -0.29  -0.59   0.34   0.31   0.60
  -0.19   0.63   0.44  -0.38   0.48
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda );

/* Parameters */
#define N 5
#define LDA N

/* Main program */
int main() {
	/* Locals */
	MKL_INT n = N, lda = LDA, info;
	/* Local arrays */
	double w[N];
	double a[LDA*N] = {
	    6.39,  0.13, -8.23, 5.71, -3.18,
	    0.00,  8.37, -4.46, -6.10,  7.21,
	    0.00,  0.00, -9.58, -9.25, -7.42,
	    0.00,  0.00, 0.00, 3.72,  8.54,
	    0.00,  0.00, 0.00, 0.00,  2.51
	};
	/* Executable statements */
	printf( "LAPACKE_dsyevd (row-major, high-level) Example Program Results\n" );
	/* Solve eigenproblem */
	info = LAPACKE_dsyevd( LAPACK_ROW_MAJOR, 'V', 'U', n, a, lda, w );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm failed to compute eigenvalues.\n" );
		exit( 1 );
	}
	/* Print eigenvalues */
	print_matrix( "Eigenvalues", 1, n, w, 1 );
	/* Print eigenvectors */
	print_matrix( "Eigenvectors (stored columnwise)", n, n, a, lda );
	exit( 0 );
} /* End of LAPACKE_dsyevd Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, double* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i*lda+j] );
		printf( "\n" );
	}
}
