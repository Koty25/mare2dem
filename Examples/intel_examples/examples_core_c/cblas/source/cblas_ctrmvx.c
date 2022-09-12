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
!      C B L A S _ C T R M V  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT          n, lda, incx;
      MKL_INT          rmaxa, cmaxa;
      MKL_Complex8    *a, *x;
      CBLAS_LAYOUT     layout;
      CBLAS_UPLO       uplo;
      CBLAS_TRANSPOSE  trans;
      CBLAS_DIAG       diag;
      MKL_INT          len_x, typeA;

      printf("\n     C B L A S _ C T R M V  EXAMPLE PROGRAM\n");

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
      if( GetIntegerParameters(in_file, &n) != 1 ) {
          printf("\n ERROR of n reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetIntegerParameters(in_file, &incx) != 1 ) {
          printf("\n ERROR of incx reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, &uplo, &trans, &diag, &layout) != 4 ) {
          printf("\n ERROR of uplo, trans, diag, layout reading\n");
          fclose(in_file);
          return 1;
      }

      rmaxa = n + 1;
      cmaxa = n;
      len_x = 1+(n-1)*ABS(incx);
      a = (MKL_Complex8 *)mkl_calloc(rmaxa*cmaxa, sizeof(MKL_Complex8), 64);
      x = (MKL_Complex8 *)mkl_calloc(len_x, sizeof(MKL_Complex8), 64);
      if( a == NULL || x == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(a);
          mkl_free(x);
          fclose(in_file);
          return 1;
      }
      if( layout == CblasRowMajor )
         lda=cmaxa;
      else
         lda=rmaxa;

      if( GetVectorC(in_file, n, x, incx) != len_x ) {
        printf("\n ERROR of vector X reading\n");
        mkl_free(a);
        mkl_free(x);
        fclose(in_file);
        return 1;
      }
      if (uplo == CblasUpper) typeA = UPPER_MATRIX;
      else                    typeA = LOWER_MATRIX;
      if( GetArrayC(in_file, &layout, typeA, &n, &n, a, &lda) != 0 ) {
        printf("\n ERROR of array A reading\n");
        mkl_free(a);
        mkl_free(x);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       N="INT_FORMAT, n);
      PrintParameters("UPLO, DIAG", uplo, diag);
      PrintParameters("TRANS", trans);
      PrintParameters("LAYOUT", layout);
      PrintVectorC(FULLPRINT, n, x, incx, "X");
      PrintArrayC(&layout, FULLPRINT, typeA, &n, &n, a, &lda, "A");

/*      Call CBLAS_CTRMV subroutine ( C Interface )            */

      cblas_ctrmv(layout, uplo, trans, diag, n, a, lda, x, incx);

/*      Print output data                                      */

      printf("\n\n     OUTPUT DATA");
      PrintVectorC(FULLPRINT, n, x, incx, "X");

      mkl_free(a);
      mkl_free(x);

      return 0;
}

