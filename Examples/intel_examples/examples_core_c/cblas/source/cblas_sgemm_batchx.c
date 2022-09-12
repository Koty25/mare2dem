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
!      C B L A S _ S G E M M _ B A T C H  Example Program Text ( C Interface )
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
      MKL_INT         m[MAX_GRP_COUNT], n[MAX_GRP_COUNT], k[MAX_GRP_COUNT];
      MKL_INT         lda[MAX_GRP_COUNT], ldb[MAX_GRP_COUNT], ldc[MAX_GRP_COUNT];
      MKL_INT         rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc;
      float           alpha[MAX_GRP_COUNT], beta[MAX_GRP_COUNT];
      float          *a, *b, *c;
      float         **a_array, **b_array, **c_array;
      CBLAS_LAYOUT    layout;
      CBLAS_TRANSPOSE transA[MAX_GRP_COUNT], transB[MAX_GRP_COUNT];
      MKL_INT         ma, na, mb, nb;
      MKL_INT         grp_count, grp_sizes[MAX_GRP_COUNT], total_batch_count, idx;

      printf("\n     C B L A S _ S G E M M _ B A T C H  EXAMPLE PROGRAM\n");

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
          if( GetIntegerParameters(in_file, &m[i], &n[i], &k[i]) != 3 ) {
              printf("\n ERROR of m, n, k reading\n");
              fclose(in_file);
              return 1;
          }
      }

      for (i=0; i<grp_count; ++i) {
          if( GetScalarsS(in_file, &alpha[i], &beta[i]) != 2 ) {
              printf("\n ERROR of alpha, beta reading\n");
              fclose(in_file);
              return 1;
          }
      }

      for (i=0; i<grp_count; ++i) {
          if( GetCblasCharParameters(in_file, &transA[i], &transB[i]) != 2 ) {
              printf("\n ERROR of transA, transB reading\n");
              fclose(in_file);
              return 1;
          }
      }

      a_array = (float **)mkl_calloc(total_batch_count, sizeof( float* ), 64);
      b_array = (float **)mkl_calloc(total_batch_count, sizeof( float* ), 64);
      c_array = (float **)mkl_calloc(total_batch_count, sizeof( float* ), 64);
      if( a_array == NULL || b_array == NULL || c_array == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(a_array);
          mkl_free(b_array);
          mkl_free(c_array);
          fclose(in_file);
          return 1;
      }

      for (l=0; l<total_batch_count; ++l) {
          a_array[l] = (float *)NULL;
          b_array[l] = (float *)NULL;
          c_array[l] = (float *)NULL;
      }

      idx = 0;
      for (i=0; i<grp_count; ++i) {
          if( transA[i] == CblasNoTrans ) {
              rmaxa = m[i] + 1;
              cmaxa = k[i];
              ma    = m[i];
              na    = k[i];
          } else {
              rmaxa = k[i] + 1;
              cmaxa = m[i];
              ma    = k[i];
              na    = m[i];
          }
          if( transB[i] == CblasNoTrans ) {
              rmaxb = k[i] + 1;
              cmaxb = n[i];
              mb    = k[i];
              nb    = n[i];
          } else {
              rmaxb = n[i] + 1;
              cmaxb = k[i];
              mb    = n[i];
              nb    = k[i];
          }
          rmaxc = m[i] + 1;
          cmaxc = n[i];

          for (j=0; j<grp_sizes[i]; ++j) {
              a = (float *)mkl_calloc(rmaxa*cmaxa, sizeof( float ), 64);
              b = (float *)mkl_calloc(rmaxb*cmaxb, sizeof( float ), 64);
              c = (float *)mkl_calloc(rmaxc*cmaxc, sizeof( float ), 64);
              if( a == NULL || b == NULL || c == NULL ) {
                  printf( "\n Can't allocate memory for arrays\n");
                  for (l=0; l<idx; ++l) {
                      mkl_free(a_array[l]);
                      mkl_free(b_array[l]);
                      mkl_free(c_array[l]);
                  }
                  mkl_free(a);
                  mkl_free(b);
                  mkl_free(c);
                  mkl_free(a_array);
                  mkl_free(b_array);
                  mkl_free(c_array);
                  fclose(in_file);
                  return 1;
              }
              if( layout == CblasRowMajor ) {
                  lda[i]=cmaxa;
                  ldb[i]=cmaxb;
                  ldc[i]=cmaxc;
              } else {
                  lda[i]=rmaxa;
                  ldb[i]=rmaxb;
                  ldc[i]=rmaxc;
              }
              if( GetArrayS(in_file, &layout, GENERAL_MATRIX, &ma, &na, a, &lda[i]) != 0 ) {
                  printf("\n ERROR of array A reading, group idx="INT_FORMAT" matrix idx="INT_FORMAT, i, idx);
                  printf("\n");
                  for (l=0; l<idx; ++l) {
                      mkl_free(a_array[l]);
                      mkl_free(b_array[l]);
                      mkl_free(c_array[l]);
                  }
                  mkl_free(a);
                  mkl_free(b);
                  mkl_free(c);
                  mkl_free(a_array);
                  mkl_free(b_array);
                  mkl_free(c_array);
                  fclose(in_file);
                  return 1;
              }
              if( GetArrayS(in_file, &layout, GENERAL_MATRIX, &mb, &nb, b, &ldb[i]) != 0 ) {
                  printf("\n ERROR of array B reading, group idx="INT_FORMAT" matrix idx="INT_FORMAT, i, idx);
                  printf("\n");
                  for (l=0; l<idx; ++l) {
                      mkl_free(a_array[l]);
                      mkl_free(b_array[l]);
                      mkl_free(c_array[l]);
                  }
                  mkl_free(a);
                  mkl_free(b);
                  mkl_free(c);
                  mkl_free(a_array);
                  mkl_free(b_array);
                  mkl_free(c_array);
                  fclose(in_file);
                  return 1;
              }
              if( GetArrayS(in_file, &layout, GENERAL_MATRIX, &m[i],  &n[i],  c, &ldc[i]) != 0 ) {
                  printf("\n ERROR of array C reading, group idx="INT_FORMAT" matrix idx="INT_FORMAT, i, idx);
                  printf("\n");
                  for (l=0; l<idx; ++l) {
                      mkl_free(a_array[l]);
                      mkl_free(b_array[l]);
                      mkl_free(c_array[l]);
                  }
                  mkl_free(a);
                  mkl_free(b);
                  mkl_free(c);
                  mkl_free(a_array);
                  mkl_free(b_array);
                  mkl_free(c_array);
                  fclose(in_file);
                  return 1;
              }

              a_array[idx] = a;
              b_array[idx] = b;
              c_array[idx] = c;

              idx++;
          }
      }
      fclose(in_file);

/*       Print input data                                      */
      printf("\n     INPUT DATA");
      PrintParameters("LAYOUT", layout);
      idx = 0;
      for (i=0; i<grp_count; ++i) {
          if( transA[i] == CblasNoTrans ) {
              ma    = m[i];
              na    = k[i];
          } else {
              ma    = k[i];
              na    = m[i];
          }
          if( transB[i] == CblasNoTrans ) {
              mb    = k[i];
              nb    = n[i];
          } else {
              mb    = n[i];
              nb    = k[i];
          }

          printf("\n       Group="INT_FORMAT", Number of matrices in the Group="INT_FORMAT, i, grp_sizes[i]);
          printf("\n       M="INT_FORMAT"  N="INT_FORMAT"  K="INT_FORMAT, m[i], n[i], k[i]);
          printf("\n       ALPHA=%5.1f  BETA=%5.1f", alpha[i], beta[i]);
          PrintParameters("TRANSA, TRANSB", transA[i], transB[i]);
          for (j=0; j<grp_sizes[i]; ++j) {
              PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &ma, &na, a_array[idx], &lda[i], "A");
              PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &mb, &nb, b_array[idx], &ldb[i], "B");
              PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &m[i],  &n[i],  c_array[idx], &ldc[i], "C");
              idx++;
          }
      }

/*      Call SGEMM_BATCH subroutine ( C Interface )                  */
      cblas_sgemm_batch(layout, transA, transB, m, n, k, alpha,
              (const float **) a_array, lda, (const float **) b_array, ldb,
              beta, c_array, ldc, grp_count, grp_sizes);

/*       Print output data                                     */
      printf("\n\n     OUTPUT DATA");
      idx = 0;
      for (i=0; i<grp_count; ++i) {
          printf("\n       Group="INT_FORMAT, i);
          for (j=0; j<grp_sizes[i]; ++j) {
              printf("\n       Matrix idx="INT_FORMAT, idx);
              PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &m[i], &n[i], c_array[idx], &ldc[i], "C");
              mkl_free(a_array[idx]);
              mkl_free(b_array[idx]);
              mkl_free(c_array[idx]);
              idx++;
          }
      }
      printf("\n");

      mkl_free(a_array);
      mkl_free(b_array);
      mkl_free(c_array);
      return 0;
}
