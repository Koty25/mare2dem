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
!      C B L A S _ G E M M _ B F 1 6 B F 1 6 F 3 2 _ C O M P U T E  
!      Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT           m, n, k;
      MKL_INT           lda, ldb, ldc;
      MKL_INT           rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc;
      float             alpha, beta;
      MKL_BF16          *a, *b, *dest;
      float             *c;
      CBLAS_LAYOUT      layout;
      CBLAS_IDENTIFIER  identifier;
      CBLAS_STORAGE     storage;
      CBLAS_TRANSPOSE   transA, transB;
      MKL_INT           ma, na, mb, nb;
      size_t            destsize;

      printf("\n     C B L A S _ G E M M _ B F 1 6 B F 1 6 F 3 2 _ C O M P U T E  EXAMPLE PROGRAM\n");

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
      if( GetIntegerParameters(in_file, &m, &n, &k) != 3 ) {
          printf("\n ERROR of m, n, k reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetScalarsS(in_file, &alpha, &beta) != 2 ) {
          printf("\n ERROR of alpha, beta reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, &transA, &transB, &layout) != 3 ) {
          printf("\n ERROR of transA, transB, layout reading\n");
          fclose(in_file);
          return 1;
      }

      if( transA == CblasNoTrans ) {
         rmaxa = m + 1;
         cmaxa = k;
         ma    = m;
         na    = k;
      } else {
         rmaxa = k + 1;
         cmaxa = m;
         ma    = k;
         na    = m;
      }
      if( transB == CblasNoTrans ) {
         rmaxb = k + 1;
         cmaxb = n;
         mb    = k;
         nb    = n;
      } else {
         rmaxb = n + 1;
         cmaxb = k;
         mb    = n;
         nb    = k;
      }
      rmaxc = m + 1;
      cmaxc = n;
      a = (MKL_BF16 *)mkl_calloc(rmaxa*cmaxa, sizeof( MKL_BF16 ), 64);
      b = (MKL_BF16 *)mkl_calloc(rmaxb*cmaxb, sizeof( MKL_BF16 ), 64);
      c = (float *)mkl_calloc(rmaxc*cmaxc, sizeof( float ), 64);
      if( a == NULL || b == NULL || c == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(a);
          mkl_free(b);
          mkl_free(c);
          fclose(in_file);
          return 1;
      }
      if (layout == CblasRowMajor) {
          lda=cmaxa;
          ldb=cmaxb;
          ldc=cmaxc;
      } else {
          lda=rmaxa;
          ldb=rmaxb;
          ldc=rmaxc;
      }
      if( GetArrayBF16(in_file, &layout, GENERAL_MATRIX, &ma, &na, a, &lda) != 0 ) {
        printf("\n ERROR of array A reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayBF16(in_file, &layout, GENERAL_MATRIX, &mb, &nb, b, &ldb) != 0 ) {
        printf("\n ERROR of array B reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayS(in_file, &layout, GENERAL_MATRIX, &m,  &n,  c, &ldc) != 0 ) {
        printf("\n ERROR of array C reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       M="INT_FORMAT"  N="INT_FORMAT"  K="INT_FORMAT, m, n, k);
      printf("\n       ALPHA=%5.1f  BETA=%5.1f", alpha, beta);
      PrintParameters("TRANSA, TRANSB", transA, transB);
      PrintParameters("LAYOUT", layout);
      PrintArrayBF16(&layout, FULLPRINT, GENERAL_MATRIX, &ma, &na, a, &lda, "A");
      PrintArrayBF16(&layout, FULLPRINT, GENERAL_MATRIX, &mb, &nb, b, &ldb, "B");
      PrintArrayS   (&layout, FULLPRINT, GENERAL_MATRIX, &m,  &n,  c, &ldc, "C");

/*      Call CBLAS_GEMM_BF16BF16F32_PACK_GET_SIZE and MKL_MALLOC functions to allocate buffer for A (C Interface) */
      identifier = CblasAMatrix;
      destsize = cblas_gemm_bf16bf16f32_pack_get_size(identifier, m, n, k);
      dest = mkl_malloc(destsize, 64);
      if( dest == NULL ) {
          printf( "\n Can't allocate memory for buffer\n");
          mkl_free(a);
          mkl_free(b);
          mkl_free(c);
          return 1;
      }
/*      Call CBLAS_GEMM_BF16BF16F32_PACK subroutine to perform scaling and packing (C Interface) */
      cblas_gemm_bf16bf16f32_pack(layout, identifier, transA, m, n, k, a, lda, dest);
/*      Call CBLAS_GEMM_BF16BF16F32_COMPUTE subroutine (C Interface)                     */
      storage = CblasPacked;
      cblas_gemm_bf16bf16f32_compute(layout, storage, transB, m, n, k,
                                     alpha, dest, lda, b, ldb, beta, c, ldc);
/*      Call MKL_FREE subroutine to free allocated storage (C Interface) */
      mkl_free(dest);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintArrayS(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, c, &ldc, "C");

      mkl_free(a);
      mkl_free(b);
      mkl_free(c);

      return 0;
}

