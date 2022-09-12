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
!      C B L A S _ D S Y R 2 K  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT         n, k;
      MKL_INT         lda, ldb, ldc;
      MKL_INT         rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc;
      double          alpha, beta;
      double         *a, *b, *c;
      CBLAS_LAYOUT    layout;
      CBLAS_TRANSPOSE trans;
      CBLAS_UPLO      uplo;
      MKL_INT         ma, na, mb, nb, typeC;

      printf("\n     C B L A S _ D S Y R 2 K  EXAMPLE PROGRAM\n");

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
      if( GetIntegerParameters(in_file, &n, &k) != 2 ) {
          printf("\n ERROR of n, k reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetScalarsD(in_file, &alpha, &beta) != 2 ) {
          printf("\n ERROR of alpha, beta reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, &uplo, &trans, &layout) != 3 ) {
          printf("\n ERROR of uplo, trans, layout reading\n");
          fclose(in_file);
          return 1;
      }

      if( trans == CblasNoTrans ) {
          rmaxa = n + 1;
          cmaxa = k;
          rmaxb = n + 1;
          cmaxb = k;
          ma    = n;
          na    = k;
          mb    = n;
          nb    = k;
      } else {
          rmaxa = k + 1;
          cmaxa = n;
          rmaxb = k + 1;
          cmaxb = n;
          ma    = k;
          na    = n;
          mb    = k;
          nb    = n;
      }
      rmaxc = n + 1;
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
      if (uplo == CblasUpper) typeC = UPPER_MATRIX;
      else                    typeC = LOWER_MATRIX;
      if( GetArrayD(in_file, &layout, typeC, &n, &n, c, &ldc) != 0 ) {
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
      printf("\n       N="INT_FORMAT"  K="INT_FORMAT, n, k);
      printf("\n       ALPHA=%6.2f  BETA=%6.2f", alpha, beta);
      PrintParameters("UPLO", uplo);
      PrintParameters("TRANS", trans);
      PrintParameters("LAYOUT", layout);
      PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &ma, &na, a, &lda, "A");
      PrintArrayD(&layout, FULLPRINT, GENERAL_MATRIX, &mb, &nb, b, &ldb, "B");
      PrintArrayD(&layout, FULLPRINT, typeC, &n, &n, c, &ldc, "C");

/*      Call CBLAS_DSYR2K subroutine ( C Interface )           */

      cblas_dsyr2k(layout, uplo, trans, n, k, alpha, a, lda,
                   b, ldb, beta, c, ldc);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintArrayD(&layout, FULLPRINT, typeC, &n, &n, c, &ldc, "C");

      mkl_free(a);
      mkl_free(b);
      mkl_free(c);

      return 0;
}

