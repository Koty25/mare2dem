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
!      C B L A S _ D T R S M _ B A T C H  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"
#include "mkl_cblas.h"

#define MAX_GRP_COUNT 5

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT         i, j, l;
      MKL_INT         m[MAX_GRP_COUNT], n[MAX_GRP_COUNT];
      MKL_INT         lda[MAX_GRP_COUNT], ldb[MAX_GRP_COUNT];
      MKL_INT         rmaxa, cmaxa, rmaxb, cmaxb;
      double          alpha[MAX_GRP_COUNT];
      double          *a, *b;
      double          **a_array, **b_array;
      CBLAS_LAYOUT    layout;
      CBLAS_SIDE      side[MAX_GRP_COUNT];
      CBLAS_UPLO      uplo[MAX_GRP_COUNT];
      CBLAS_TRANSPOSE transA[MAX_GRP_COUNT];
      CBLAS_DIAG      diag[MAX_GRP_COUNT];
      MKL_INT         ma, na, typeA;
      MKL_INT         grp_count, grp_sizes[MAX_GRP_COUNT], total_batch_count, idx;

      printf("\n     C B L A S _ D T R S M _ B A T C H  EXAMPLE PROGRAM\n");

/*       Get input parameters                                  */

      if( argc == 1 ) {
          printf("\n You must specify in_file data file as 1-st parameter");
          return 1;
      }
      in_file_name = argv[1];

/*       Get input data                                        */

      if( (in_file = fopen( in_file_name, "r" )) == NULL ) {
         printf("\n ERROR on OPEN '%s' with mode=\"r\"\n", in_file_name);
         return 1;
      }

      if( GetCblasCharParameters(in_file, &layout) != 1 ) {
          printf("\n ERROR of layout reading\n");
          fclose(in_file);
          return 1;
      }

      if( GetIntegerParameters(in_file, &grp_count) != 1 ) {
          printf("\n ERROR of grp_count reading\n");
          fclose(in_file);
          return 1;
      }

      if ( grp_count > MAX_GRP_COUNT ) {
          printf("\n grp_count exceeds the limit\n");
          fclose(in_file);
          return 1;
      }

      if( GetIntegerParameters(in_file, &grp_sizes[0], &grp_sizes[1], &grp_sizes[2],
                  &grp_sizes[3], &grp_sizes[4]) != grp_count ) {
          printf("\n ERROR of grp_sizes reading\n");
          fclose(in_file);
          return 1;
      }

      total_batch_count = 0;
      for (i=0; i<grp_count; ++i) {
          total_batch_count += grp_sizes[i];
          if( GetIntegerParameters(in_file, &m[i], &n[i]) != 2 ) {
              printf("\n ERROR of m, n reading\n");
              fclose(in_file);
              return 1;
          }
      }

      for (i=0; i<grp_count; ++i) {
          if( GetScalarsD(in_file, &alpha[i]) != 1 ) {
              printf("\n ERROR of alpha, beta reading\n");
              fclose(in_file);
              return 1;
          }
      }

      for (i=0; i<grp_count; ++i) {
          if( GetCblasCharParameters(in_file, &side[i], &uplo[i], &transA[i], &diag[i]) != 4 ) {
              printf("\n ERROR of side, uplo, transA, diag reading\n");
              fclose(in_file);
              return 1;
          }
      }

      a_array = (double **)mkl_calloc(total_batch_count, sizeof( double* ), 64);
      b_array = (double **)mkl_calloc(total_batch_count, sizeof( double* ), 64);
      if( a_array == NULL || b_array == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(a_array);
          mkl_free(b_array);
          fclose(in_file);
          return 1;
      }

      for (l=0; l<total_batch_count; ++l) {
          a_array[l] = (double *)NULL;
          b_array[l] = (double *)NULL;
      }

      idx = 0;
      for (i=0; i<grp_count; ++i) {
          if( side[i] == CblasLeft ) {
              rmaxa = m[i] + 1;
              cmaxa = m[i];
              ma    = m[i];
              na    = m[i];
          } else {
              rmaxa = n[i] + 1;
              cmaxa = n[i];
              ma    = n[i];
              na    = n[i];
          }
          rmaxb = m[i] + 1;
          cmaxb = n[i];
          if (uplo[i] == CblasUpper) typeA = UPPER_MATRIX;
          else                       typeA = LOWER_MATRIX;

          if( layout == CblasRowMajor ) {
              lda[i]=cmaxa;
              ldb[i]=cmaxb;
          } else {
              lda[i]=rmaxa;
              ldb[i]=rmaxb;
          }

          for (j=0; j<grp_sizes[i]; ++j) {
              a = (double *)mkl_calloc(rmaxa*cmaxa, sizeof( double ), 64);
              b = (double *)mkl_calloc(rmaxb*cmaxb, sizeof( double ), 64);
              if( a == NULL || b == NULL ) {
                  printf( "\n Can't allocate memory for arrays\n");
                  for (l=0; l<idx; ++l) {
                      mkl_free(a_array[l]);
                      mkl_free(b_array[l]);
                  }
                  mkl_free(a);
                  mkl_free(b);
                  mkl_free(a_array);
                  mkl_free(b_array);
                  fclose(in_file);
                  return 1;
              }

              if( GetArrayD(in_file, &layout, typeA, &ma, &na, a, &lda[i]) != 0 ) {
                  printf("\n ERROR of array A reading, group idx="INT_FORMAT" matrix idx="INT_FORMAT, i, idx);
                  printf("\n");
                  for (l=0; l<idx; ++l) {
                      mkl_free(a_array[l]);
                      mkl_free(b_array[l]);
                  }
                  mkl_free(a);
                  mkl_free(b);
                  mkl_free(a_array);
                  mkl_free(b_array);
                  fclose(in_file);
                  return 1;
              }
              if( GetArrayD(in_file, &layout, GENERAL_MATRIX, &m[i], &n[i], b, &ldb[i]) != 0 ) {
                  printf("\n ERROR of array B reading, group idx="INT_FORMAT" matrix idx="INT_FORMAT, i, idx);
                  printf("\n");
                  for (l=0; l<idx; ++l) {
                      mkl_free(a_array[l]);
                      mkl_free(b_array[l]);
                  }
                  mkl_free(a);
                  mkl_free(b);
                  mkl_free(a_array);
                  mkl_free(b_array);
                  fclose(in_file);
                  return 1;
              }

              a_array[idx] = a;
              b_array[idx] = b;

              idx++;
          }
      }
      fclose(in_file);

/*       Print input data                                      */
      printf("\n     INPUT DATA");
      PrintParameters("LAYOUT", layout);
      idx = 0;
      for (i=0; i<grp_count; ++i) {
          if( side[i] == CblasLeft ) {
              ma    = m[i];
              na    = m[i];
          } else {
              ma    = n[i];
              na    = n[i];
          }
          if (uplo[i] == CblasUpper) typeA = UPPER_MATRIX;
          else                       typeA = LOWER_MATRIX;
          printf("\n       Group="INT_FORMAT", Number of matrices in the Group="INT_FORMAT, i, grp_sizes[i]);
          printf("\n       M="INT_FORMAT"  N="INT_FORMAT, m[i], n[i]);
          printf("\n       ALPHA=%5.1f", alpha[i]);
          PrintParameters("SIDE, UPLO, TRANSA, DIAG", side[i], uplo[i], transA[i], diag[i]);
          for (j=0; j<grp_sizes[i]; ++j) {
              PrintArrayD(&layout, FULLPRINT, typeA, &ma, &na, a_array[idx], &lda[i], "A");
              PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &m[i], &n[i], b_array[idx], &ldb[i], "B");
              idx++;
          }
      }

/*      Call DTRSM_BATCH subroutine ( C Interface )                  */
      cblas_dtrsm_batch(layout, side, uplo, transA, diag, m, n, alpha,
              (const double **) a_array, lda, b_array, ldb, grp_count,
              grp_sizes);

/*       Print output data                                     */
      printf("\n\n     OUTPUT DATA");
      idx = 0;
      for (i=0; i<grp_count; ++i) {
          printf("\n       Group="INT_FORMAT, i);
          for (j=0; j<grp_sizes[i]; ++j) {
              printf("\n       Matrix idx="INT_FORMAT, idx);
              PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &m[i], &n[i], b_array[idx], &ldb[i], "B");
              mkl_free(a_array[idx]);
              mkl_free(b_array[idx]);
              idx++;
          }
      }
      printf("\n");

      mkl_free(a_array);
      mkl_free(b_array);
      return 0;
}
