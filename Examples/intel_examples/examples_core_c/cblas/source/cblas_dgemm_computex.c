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
!      C B L A S _ D G E M M _ C O M P U T E  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"
#include "mkl_cblas.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT           m, n, k;
      MKL_INT           lda, ldb, ldc;
      MKL_INT           rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc;
      double            alpha, beta;
      double            *a, *b, *c, *dest;
      CBLAS_LAYOUT      layout;
      CBLAS_IDENTIFIER  identifier;
      CBLAS_STORAGE     storage;
      CBLAS_TRANSPOSE   transA, transB;
      MKL_INT           ma, na, mb, nb;
      size_t            destsize;

      printf("\n     C B L A S _ D G E M M _ C O M P U T E  EXAMPLE PROGRAM\n");

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
      if( GetScalarsD(in_file, &alpha, &beta) != 2 ) {
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
      a = (double *)mkl_calloc(rmaxa*cmaxa, sizeof( double ), 64);
      b = (double *)mkl_calloc(rmaxb*cmaxb, sizeof( double ), 64);
      c = (double *)mkl_calloc(rmaxc*cmaxc, sizeof( double ), 64);
      if( a == NULL || b == NULL || c == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(a);
          mkl_free(b);
          mkl_free(c);
          fclose(in_file);
          return 1;
      }
      if( layout == CblasRowMajor ) {
          lda=cmaxa;
          ldb=cmaxb;
          ldc=cmaxc;
      } else {
          lda=rmaxa;
          ldb=rmaxb;
          ldc=rmaxc;
      }
      if( GetArrayD(in_file, &layout, GENERAL_MATRIX, &ma, &na, a, &lda) != 0 ) {
        printf("\n ERROR of array A reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayD(in_file, &layout, GENERAL_MATRIX, &mb, &nb, b, &ldb) != 0 ) {
        printf("\n ERROR of array B reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayD(in_file, &layout, GENERAL_MATRIX, &m,  &n,  c, &ldc) != 0 ) {
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
      PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &ma, &na, a, &lda, "A");
      PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &mb, &nb, b, &ldb, "B");
      PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &m,  &n,  c, &ldc, "C");

/*      Call CBLAS_DGEMM_PACK_GET_SIZE and MKL_MALLOC functions to allocate buffer for B (C Interface) */
      identifier = CblasBMatrix;
      destsize = cblas_dgemm_pack_get_size(identifier, m, n, k);
      dest = mkl_malloc(destsize, 64);
      if( dest == NULL ) {
          printf( "\n Can't allocate memory for buffer\n");
          mkl_free(a);
          mkl_free(b);
          mkl_free(c);
          return 1;
      }
/*      Call CBLAS_DGEMM_PACK subroutine to perform scaling and packing (C Interface) */
      cblas_dgemm_pack(layout, identifier, transB, m, n, k, alpha, b, ldb, dest);
/*      Call CBLAS_DGEMM_COMPUTE subroutine (C Interface)                     */
      storage = CblasPacked;
      cblas_dgemm_compute(layout, transA, storage, m, n, k,
                  a, lda, dest, ldb, beta, c, ldc);
/*      Call MKL_FREE subroutine to free allocated storage (C Interface) */
      mkl_free(dest);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, c, &ldc, "C");

      mkl_free(a);
      mkl_free(b);
      mkl_free(c);

      return 0;
}

