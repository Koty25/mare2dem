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
!      C B L A S _ Z H E R 2 K  Example Program Text ( C Interface )
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
      MKL_Complex16   alpha;
      double          beta;
      MKL_Complex16  *a, *b, *c;
      CBLAS_LAYOUT    layout;
      CBLAS_TRANSPOSE trans;
      CBLAS_UPLO      uplo;
      MKL_INT         ma, na, mb, nb, typeC;

      printf("\n     C B L A S _ Z H E R 2 K  EXAMPLE PROGRAM\n");

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
      if( GetScalarsZ(in_file, &alpha) != 1 ) {
          printf("\n ERROR of alpha reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetScalarsD(in_file, &beta) != 1 ) {
          printf("\n ERROR of beta reading\n");
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
      a = (MKL_Complex16 *)mkl_calloc(rmaxa*cmaxa, sizeof(MKL_Complex16), 64);
      b = (MKL_Complex16 *)mkl_calloc(rmaxb*cmaxb, sizeof(MKL_Complex16), 64);
      c = (MKL_Complex16 *)mkl_calloc(rmaxc*cmaxc, sizeof(MKL_Complex16), 64);
      if ( a == NULL || b == NULL || c == NULL ) {
          printf("\n Can't allocate memory arrays");
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
      if( GetArrayZ(in_file, &layout, 0, &ma, &na, a, &lda) != 0 ) {
        printf("\n ERROR of array A reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayZ(in_file, &layout, 0, &mb, &nb, b, &ldb) != 0 ) {
        printf("\n ERROR of array B reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if (uplo == CblasUpper) typeC = UPPER_MATRIX;
      else                    typeC = LOWER_MATRIX;
      if( GetArrayZ(in_file, &layout, typeC, &n, &n, c, &ldc) != 0 ) {
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
      printf("\n       ALPHA =(%4.1f,%4.1f)  BETA =%4.1f",
              alpha.real, alpha.imag, beta);
      PrintParameters("UPLO", uplo);
      PrintParameters("TRANS", trans);
      PrintParameters("LAYOUT", layout);
      PrintArrayZ(&layout, FULLPRINT, GENERAL_MATRIX, &ma, &na, a, &lda, "A");
      PrintArrayZ(&layout, FULLPRINT, GENERAL_MATRIX, &mb, &nb, b, &ldb, "B");
      PrintArrayZ(&layout, FULLPRINT, typeC, &n, &n, c, &ldc, "C");

/*      Call CBLAS_ZHER2K subroutine ( C Interface )           */

      cblas_zher2k(layout, uplo, trans, n, k, &alpha, a, lda,
                   b, ldb, beta, c, ldc);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintArrayZ(&layout, FULLPRINT, typeC, &n, &n, c, &ldc, "C");

      mkl_free(a);
      mkl_free(b);
      mkl_free(c);

      return 0;
}

