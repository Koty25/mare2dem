/*******************************************************************************
* Copyright 1999-2020 Intel Corporation.
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
!  Content:
!      M K L _ D R T S M _ C O M P A C T  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

#define MEM_ALIGNMENT 64
#define MAX_NM 512

int main(int argc, char *argv[])
{

  FILE *in_file;
  char *in_file_name;

  MKL_INT i, j, p_num, idx_matrix;
  MKL_INT m, n, nm;
  MKL_INT lda, ldap, ldb, ldbp;
  MKL_INT sdap, sdbp;
  MKL_INT a_idx, b_idx, adim;
  MKL_INT a_total, b_total;
  double *a, *b;
  double *a_array[MAX_NM], *b_array[MAX_NM];
  double *a_compact, *b_compact;
  MKL_INT a_buffer, b_buffer;
  CBLAS_SIDE SIDE;
  CBLAS_UPLO UPLO;
  CBLAS_DIAG DIAG;
  CBLAS_TRANSPOSE TRANSA;
  CBLAS_LAYOUT layout;
  MKL_COMPACT_PACK COMPACT_FORMAT;
  double alpha;

  printf("\n     M K L _ D T R S M _ C O M P A C T  EXAMPLE PROGRAM\n");

/*       Get input data                                        */

  if( argc == 1 ) {
    printf("\n You must specify in_file data file as 1-st parameter");
    return 1;
  }
  in_file_name = argv[1];

  if( (in_file = fopen( in_file_name, "r" )) == NULL ) {
     printf("\n ERROR on OPEN '%s' with mode=\"r\"\n", in_file_name);
     return 1;
  }
  if( GetIntegerParameters(in_file, &m, &n, &nm) != 3 ) {
      printf("\n ERROR of m, n, k reading\n");
      fclose(in_file);
      return 1;
  }
  if( GetScalarsD(in_file, &alpha) != 1 ) {
      printf("\n ERROR of alpha, beta reading\n");
      fclose(in_file);
      return 1;
  }
  if( GetCblasCharParameters(in_file, &SIDE, &UPLO, &TRANSA, &DIAG, &layout) !=5 ) {
    printf("\n ERROR of SIDE, UPLO, TRANSA, DIAG, layout reading\n");
    fclose(in_file);
    return 1;
  }
  fclose(in_file);

  if ( SIDE == CblasLeft ) {
    adim = m;
  } else {
    adim = n;
  }
  lda = adim;
  ldap = adim;
  sdap = adim;

  if ( layout == CblasColMajor ) {
    ldb = m;
    ldbp = m;
    sdbp = n;
  } else {
    ldb = n;
    ldbp = n;
    sdbp = m;
  }

/*	Set up standard arrays in P2P format			*/
  a_total = adim * adim * nm;
  b_total = m * n * nm;
  a = (double *)mkl_malloc( a_total * sizeof(double), MEM_ALIGNMENT );
  b = (double *)mkl_malloc( b_total * sizeof(double), MEM_ALIGNMENT );

  for (i = 0; i < a_total; i++) a[i] = rand() / (double) RAND_MAX + .5;
  for (i = 0; i < b_total; i++) b[i] = rand() / (double) RAND_MAX + .5;

  a_idx = 0;
  b_idx = 0;
  p_num = 0;
  for (j = 0; j < nm; j++) {
    a_array[p_num] = &a[ a_idx ];
    b_array[p_num] = &b[ b_idx ];
    p_num++;
    a_idx += adim * lda;
    b_idx += (layout == CblasColMajor) ? (ldb * n) : (ldb * m);
  }
  idx_matrix = 0;
  for (j = 0; j < nm; j++) {
    for (i = 0; i < adim; i++) {
      a_array[idx_matrix][i*lda + i] = 10.0;
    }
    idx_matrix++;
  }

/*       Print input data                                      */

  printf("\n     INPUT DATA");
  printf("\n       M="INT_FORMAT"  N="INT_FORMAT"  NM="INT_FORMAT, m, n, nm);
  printf("\n       ALPHA=%5.1f", alpha);
  PrintParameters("SIDE, UPLO, TRANSA, DIAG", SIDE, UPLO, TRANSA, DIAG);
  PrintParameters("LAYOUT", layout);
  for (i = 0; i < nm; i++) {
    printf("\nMatrix set "INT_FORMAT":", i);
    PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &adim, &adim, a_array[i], &lda, "A");
    PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, b_array[i], &ldb, "B");
  }

/*	Set up Compact arrays					*/
  COMPACT_FORMAT = mkl_get_format_compact();

  a_buffer = mkl_dget_size_compact( ldap, sdap, COMPACT_FORMAT, nm );
  b_buffer = mkl_dget_size_compact( ldbp, sdbp, COMPACT_FORMAT, nm );

  a_compact = (double *)mkl_malloc( a_buffer, MEM_ALIGNMENT );
  b_compact = (double *)mkl_malloc( b_buffer, MEM_ALIGNMENT );

/*	Pack from P2P to Compact format				*/
  mkl_dgepack_compact( (MKL_LAYOUT)layout, adim, adim, (const double* const*)a_array, lda, a_compact, ldap, COMPACT_FORMAT, nm );
  mkl_dgepack_compact( (MKL_LAYOUT)layout, m, n, (const double* const*)b_array, ldb, b_compact, ldbp, COMPACT_FORMAT, nm );

/*	Perform Compact GEMM					*/
  mkl_dtrsm_compact( (MKL_LAYOUT)layout, (MKL_SIDE)SIDE, (MKL_UPLO)UPLO, (MKL_TRANSPOSE)TRANSA, (MKL_DIAG)DIAG, m, n, alpha, a_compact, ldap, b_compact, ldbp, COMPACT_FORMAT, nm );

/*	Unpack from Compact to P2P format			*/
  mkl_dgeunpack_compact( (MKL_LAYOUT)layout, m, n, b_array, ldb, b_compact, ldbp, COMPACT_FORMAT, nm );

/*       Print output data                                     */
  printf("\n\n     OUTPUT DATA");
  for (i = 0; i < nm; i++) {
    printf("\nMatrix set "INT_FORMAT":", i);
    PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &m,  &n,  b_array[i], &ldb, "B");
  }

  mkl_free(a_compact);
  mkl_free(b_compact);
  mkl_free(a);
  mkl_free(b);

  return 0;
}
