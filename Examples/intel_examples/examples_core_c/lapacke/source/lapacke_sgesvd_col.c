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
   LAPACKE_sgesvd Example.
   =======================

   Program computes the singular value decomposition of a general
   rectangular matrix A:

     8.79   9.93   9.83   5.45   3.16
     6.11   6.91   5.04  -0.27   7.98
    -9.15  -7.93   4.86   4.85   3.01
     9.57   1.64   8.83   0.74   5.80
    -3.49   4.02   9.80  10.00   4.27
     9.84   0.15  -8.99  -6.02  -5.31

   Description.
   ============

   The routine computes the singular value decomposition (SVD) of a real
   m-by-n matrix A, optionally computing the left and/or right singular
   vectors. The SVD is written as

   A = U*SIGMA*VT

   where SIGMA is an m-by-n matrix which is zero except for its min(m,n)
   diagonal elements, U is an m-by-m orthogonal matrix and VT (V transposed)
   is an n-by-n orthogonal matrix. The diagonal elements of SIGMA
   are the singular values of A; they are real and non-negative, and are
   returned in descending order. The first min(m, n) columns of U and V are
   the left and right singular vectors of A.

   Note that the routine returns VT, not V.

   Example Program Results.
   ========================

 LAPACKE_sgesvd (column-major, high-level) Example Program Results

 Singular values
  27.47  22.64   8.56   5.99   2.01

 Left singular vectors (stored columnwise)
  -0.59   0.26   0.36   0.31   0.23
  -0.40   0.24  -0.22  -0.75  -0.36
  -0.03  -0.60  -0.45   0.23  -0.31
  -0.43   0.24  -0.69   0.33   0.16
  -0.47  -0.35   0.39   0.16  -0.52
   0.29   0.58  -0.02   0.38  -0.65

 Right singular vectors (stored rowwise)
  -0.25  -0.40  -0.69  -0.37  -0.41
   0.81   0.36  -0.25  -0.37  -0.10
  -0.26   0.70  -0.22   0.39  -0.49
   0.40  -0.45   0.25   0.43  -0.62
  -0.22   0.14   0.59  -0.63  -0.44
*/
#include <stdlib.h>
#include <stdio.h>
#include "mkl_lapacke.h"

#define min(a,b) ((a)>(b)?(b):(a))

/* Auxiliary routines prototypes */
extern void print_matrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda );

/* Parameters */
#define M 6
#define N 5
#define LDA M
#define LDU M
#define LDVT N

/* Main program */
int main() {
	/* Locals */
	MKL_INT m = M, n = N, lda = LDA, ldu = LDU, ldvt = LDVT, info;
	float superb[min(M,N)-1];
	/* Local arrays */
	float s[N], u[LDU*M], vt[LDVT*N];
	float a[LDA*N] = {
	    8.79f,  6.11f, -9.15f,  9.57f, -3.49f,  9.84f,
	    9.93f,  6.91f, -7.93f,  1.64f,  4.02f,  0.15f,
	    9.83f,  5.04f,  4.86f,  8.83f,  9.80f, -8.99f,
	    5.45f, -0.27f,  4.85f,  0.74f, 10.00f, -6.02f,
	    3.16f,  7.98f,  3.01f,  5.80f,  4.27f, -5.31f
	};
	/* Executable statements */
	printf( "LAPACKE_sgesvd (column-major, high-level) Example Program Results\n" );
	/* Compute SVD */
	info = LAPACKE_sgesvd( LAPACK_COL_MAJOR, 'A', 'A', m, n, a, lda,
			s, u, ldu, vt, ldvt, superb );
	/* Check for convergence */
	if( info > 0 ) {
		printf( "The algorithm computing SVD failed to converge.\n" );
		exit( 1 );
	}
	/* Print singular values */
	print_matrix( "Singular values", 1, n, s, 1 );
	/* Print left singular vectors */
	print_matrix( "Left singular vectors (stored columnwise)", m, n, u, ldu );
	/* Print right singular vectors */
	print_matrix( "Right singular vectors (stored rowwise)", n, n, vt, ldvt );
	exit( 0 );
} /* End of LAPACKE_sgesvd Example */

/* Auxiliary routine: printing a matrix */
void print_matrix( char* desc, MKL_INT m, MKL_INT n, float* a, MKL_INT lda ) {
	MKL_INT i, j;
	printf( "\n %s\n", desc );
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < n; j++ ) printf( " %6.2f", a[i+j*lda] );
		printf( "\n" );
	}
}
