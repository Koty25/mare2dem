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
!      M K L _ S G E M M _ C O M P A C T  Example Program Text ( C Interface )
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

  MKL_INT i, j, p_num;
  MKL_INT m, n, k, nm;
  MKL_INT lda, ldap, ldb, ldbp, ldc, ldcp;
  MKL_INT sdap, sdbp, sdcp;
  MKL_INT a_idx, b_idx, c_idx, arows, acols, brows, bcols;
  MKL_INT a_total, b_total, c_total;
  float *a, *b, *c;
  float *a_array[MAX_NM], *b_array[MAX_NM], *c_array[MAX_NM];
  float *a_compact, *b_compact, *c_compact;
  MKL_INT a_buffer, b_buffer, c_buffer;
  CBLAS_TRANSPOSE TRANSA, TRANSB;
  CBLAS_LAYOUT layout;
  MKL_COMPACT_PACK COMPACT_FORMAT;
  float alpha, beta;

  printf("\n     M K L _ S G E M M _ C O M P A C T  EXAMPLE PROGRAM\n");

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
  if( GetIntegerParameters(in_file, &m, &n, &k, &nm) != 4 ) {
      printf("\n ERROR of m, n, k reading\n");
      fclose(in_file);
      return 1;
  }
  if( GetScalarsS(in_file, &alpha, &beta) != 2 ) {
      printf("\n ERROR of alpha, beta reading\n");
      fclose(in_file);
      return 1;
  }
  if( GetCblasCharParameters(in_file, &TRANSA, &TRANSB, &layout) !=3 ) {
    printf("\n ERROR of TRANSA, TRANSB, layout reading\n");
    fclose(in_file);
    return 1;
  }
  fclose(in_file);

  if ( TRANSA == CblasNoTrans ) {
    arows = m;
    acols = k;
  } else {
    arows = k;
    acols = m;
  }
  if ( TRANSB == CblasNoTrans ) {
    brows = k;
    bcols = n;
  } else {
    brows = n;
    bcols = k;
  }

  if ( layout == CblasColMajor ) {
    lda = arows;
    ldap = arows;
    sdap = acols;
    ldb = brows;
    ldbp = brows;
    sdbp = bcols;
    ldc = m;
    ldcp = m;
    sdcp = n;
  } else {
    lda = acols;
    ldap = acols;
    sdap = arows;
    ldb = bcols;
    ldbp = bcols;
    sdbp = brows;
    ldc = n;
    ldcp = n;
    sdcp = m;
  }

/*	Set up standard arrays in P2P format			*/
  a_total = arows * acols * nm;
  b_total = brows * bcols * nm;
  c_total = m * n * nm;
  a = (float *)mkl_malloc( a_total * sizeof(float), MEM_ALIGNMENT );
  b = (float *)mkl_malloc( b_total * sizeof(float), MEM_ALIGNMENT );
  c = (float *)mkl_malloc( c_total * sizeof(float), MEM_ALIGNMENT );

  for (i = 0; i < a_total; i++) a[i] = rand() / (float) RAND_MAX + .5;
  for (i = 0; i < b_total; i++) b[i] = rand() / (float) RAND_MAX + .5;
  for (i = 0; i < c_total; i++) c[i] = rand() / (float) RAND_MAX + .5;

  a_idx = 0;
  b_idx = 0;
  c_idx = 0;
  p_num = 0;
  for (j = 0; j < nm; j++) {
    a_array[p_num] = &a[ a_idx ];
    b_array[p_num] = &b[ b_idx ];
    c_array[p_num] = &c[ c_idx ];
    p_num++;
    a_idx += (layout == CblasColMajor) ? (lda * acols) : (lda * arows);
    b_idx += (layout == CblasColMajor) ? (ldb * bcols) : (ldb * brows);
    c_idx += (layout == CblasColMajor) ? (ldc * n) : (ldc * m);
  }

/*       Print input data                                      */

  printf("\n     INPUT DATA");
  printf("\n       M="INT_FORMAT"  N="INT_FORMAT"  K="INT_FORMAT"  NM="INT_FORMAT, m, n, k, nm);
  printf("\n       ALPHA=%5.1f  BETA=%5.1f", alpha, beta);
  PrintParameters("TRANSA, TRANSB", TRANSA, TRANSB);
  PrintParameters("LAYOUT", layout);
  for (i = 0; i < nm; i++) {
    printf("\nMatrix set "INT_FORMAT":", i);
    PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &arows, &acols, a_array[i], &lda, "A");
    PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &brows, &bcols, b_array[i], &ldb, "B");
    PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &m,  &n,  c_array[i], &ldc, "C");
  }

/*	Set up Compact arrays					*/
  COMPACT_FORMAT = mkl_get_format_compact();

  a_buffer = mkl_sget_size_compact( ldap, sdap, COMPACT_FORMAT, nm );
  b_buffer = mkl_sget_size_compact( ldbp, sdbp, COMPACT_FORMAT, nm );
  c_buffer = mkl_sget_size_compact( ldcp, sdcp, COMPACT_FORMAT, nm );

  a_compact = (float *)mkl_malloc( a_buffer, MEM_ALIGNMENT );
  b_compact = (float *)mkl_malloc( b_buffer, MEM_ALIGNMENT );
  c_compact = (float *)mkl_malloc( c_buffer, MEM_ALIGNMENT );

/*	Pack from P2P to Compact format				*/
  mkl_sgepack_compact( (MKL_LAYOUT)layout, arows, acols, (const float* const*)a_array, lda, a_compact, ldap, COMPACT_FORMAT, nm );
  mkl_sgepack_compact( (MKL_LAYOUT)layout, brows, bcols, (const float* const*)b_array, ldb, b_compact, ldbp, COMPACT_FORMAT, nm );
  mkl_sgepack_compact( (MKL_LAYOUT)layout, m, n, (const float* const*)c_array, ldc, c_compact, ldcp, COMPACT_FORMAT, nm );

/*	Perform Compact GEMM					*/
  mkl_sgemm_compact( (MKL_LAYOUT)layout, (MKL_TRANSPOSE)TRANSA, (MKL_TRANSPOSE)TRANSB, m, n, k, alpha, a_compact, ldap, b_compact, ldbp, beta, c_compact, ldcp, COMPACT_FORMAT, nm );

/*	Unpack from Compact to P2P format			*/
  mkl_sgeunpack_compact( (MKL_LAYOUT)layout, m, n, c_array, ldc, c_compact, ldcp, COMPACT_FORMAT, nm );

/*       Print output data                                     */
  printf("\n\n     OUTPUT DATA");
  for (i = 0; i < nm; i++) {
    printf("\nMatrix set "INT_FORMAT":", i);
    PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &m,  &n,  c_array[i], &ldc, "C");
  }

  mkl_free(a_compact);
  mkl_free(b_compact);
  mkl_free(c_compact);
  mkl_free(a);
  mkl_free(b);
  mkl_free(c);

  return 0;
}
