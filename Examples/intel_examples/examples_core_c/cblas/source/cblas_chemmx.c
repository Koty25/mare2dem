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
!      C B L A S _ C H E M M  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT         m, n;
      MKL_INT         lda, ldb, ldc;
      MKL_INT         rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc;
      MKL_Complex8    alpha, beta;
      MKL_Complex8   *a, *b, *c;
      CBLAS_LAYOUT    layout;
      CBLAS_SIDE      side;
      CBLAS_UPLO      uplo;
      MKL_INT         ma, na, typeA;

      printf("\n     C B L A S _ C H E M M  EXAMPLE PROGRAM\n");

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
      if( GetIntegerParameters(in_file, &m, &n) != 2 ) {
          printf("\n ERROR of m, n reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetScalarsC(in_file, &alpha, &beta) != 2 ) {
          printf("\n ERROR of alpha, beta reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, &side, &uplo, &layout) != 3 ) {
          printf("\n ERROR of side, uplo, layout reading\n");
          fclose(in_file);
          return 1;
      }

      if( side == CblasLeft ) {
          rmaxa = m + 1;
          cmaxa = m;
          ma    = m;
          na    = m;
      } else {
          rmaxa = n + 1;
          cmaxa = n;
          ma    = n;
          na    = n;
      }
      rmaxb = m + 1;
      cmaxb = n;
      rmaxc = m + 1;
      cmaxc = n;
      a = (MKL_Complex8 *)mkl_calloc(rmaxa*cmaxa, sizeof( MKL_Complex8 ), 64);
      b = (MKL_Complex8 *)mkl_calloc(rmaxb*cmaxb, sizeof( MKL_Complex8 ), 64);
      c = (MKL_Complex8 *)mkl_calloc(rmaxc*cmaxc, sizeof( MKL_Complex8 ), 64);
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
      if( uplo == CblasUpper ) typeA = UPPER_MATRIX;
      else                     typeA = LOWER_MATRIX;

      if( GetArrayC(in_file, &layout, typeA, &ma, &na, a, &lda) != 0 ) {
        printf("\n ERROR of array A reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayC(in_file, &layout, GENERAL_MATRIX, &m, &n, b, &ldb) != 0 ) {
        printf("\n ERROR of array B reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayC(in_file, &layout, GENERAL_MATRIX, &m, &n, c, &ldc) != 0 ) {
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
      printf("\n       M="INT_FORMAT"  N="INT_FORMAT, m, n);
      printf("\n       ALPHA =(%4.1f,%4.1f)  BETA =(%4.1f,%4.1f)",
              alpha.real, alpha.imag, beta.real, beta.imag);
      PrintParameters("SIDE, UPLO", side, uplo);
      PrintParameters("LAYOUT", layout);
      PrintArrayC(&layout, FULLPRINT, typeA, &ma, &na, a, &lda, "A");
      PrintArrayC(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, b, &ldb, "B");
      PrintArrayC(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, c, &ldc, "C");

/*      Call CBLAS_CHEMM subroutine ( C Interface )            */

      cblas_chemm(layout, side, uplo, m, n, &alpha, a, lda,
                  b, ldb, &beta, c, ldc);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintArrayC(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, c, &ldc, "C");

      mkl_free(a);
      mkl_free(b);
      mkl_free(c);

      return 0;
}

