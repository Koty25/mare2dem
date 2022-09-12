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
!      C B L A S _ Z H B M V  Example Program Text ( C Interface )
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

      MKL_INT         n, k, lda, incx, incy;
      MKL_INT         rmaxa, cmaxa;
      MKL_Complex16   alpha, beta;
      MKL_Complex16  *a, *x, *y;
      CBLAS_LAYOUT    layout;
      CBLAS_UPLO      uplo;
      MKL_INT         kl, ku;
      MKL_INT         len_x, len_y;

      printf("\n     C B L A S _ Z H B M V  EXAMPLE PROGRAM\n");

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
      if( GetIntegerParameters(in_file, &incx, &incy) != 2 ) {
          printf("\n ERROR of incx, incy reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetScalarsZ(in_file, &alpha, &beta) != 2 ) {
          printf("\n ERROR of alpha, beta reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, &uplo, &layout) != 2 ) {
          printf("\n ERROR of uplo, layout reading\n");
          fclose(in_file);
          return 1;
      }

      rmaxa = k + 1;
      cmaxa = n;
      len_x = 1+(n-1)*ABS(incx);
      len_y = 1+(n-1)*ABS(incy);
      a = (MKL_Complex16 *)mkl_calloc(rmaxa*cmaxa, sizeof(MKL_Complex16), 64);
      x = (MKL_Complex16 *)mkl_calloc(len_x, sizeof(MKL_Complex16), 64);
      y = (MKL_Complex16 *)mkl_calloc(len_y, sizeof(MKL_Complex16), 64);
      if( a == NULL || x == NULL || y == NULL ) {
          printf("\n Can't allocate memory for arrays\n");
          mkl_free(a);
          mkl_free(x);
          mkl_free(y);
          fclose(in_file);
          return 1;
      }
      if( GetVectorZ(in_file, n, x, incx) != len_x ) {
        printf("\n ERROR of vector X reading\n");
        mkl_free(a);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      if( GetVectorZ(in_file, n, y, incy) != len_y ) {
        printf("\n ERROR of vector Y reading\n");
        mkl_free(a);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      lda=rmaxa;
      if( uplo == CblasUpper ) {
        kl = 0;
        ku = k;
      } else {
        kl = k;
        ku = 0;
      }
      if( GetBandArrayZ(in_file, &layout, kl, ku, n, n, a, lda) != 0 ) {
        printf("\n ERROR of array A reading\n");
        mkl_free(a);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       N="INT_FORMAT"  K="INT_FORMAT, n, k);
      printf("\n       ALPHA =(%4.1f, %4.1f )  BETA =(%4.1f, %4.1f )",
              alpha.real, alpha.imag, beta.real, beta.imag);
      PrintParameters("UPLO", uplo);
      PrintParameters("LAYOUT", layout);
      PrintVectorZ(FULLPRINT, n, x, incx, "X");
      PrintVectorZ(FULLPRINT, n, y, incy, "Y");
      PrintBandArrayZ(&layout, SHORTPRINT, kl, ku, n, n, a, lda, "A");

/*      Call CBLAS_ZHBMV subroutine ( C Interface )            */

      cblas_zhbmv(layout, uplo, n, k, &alpha, a, lda, x, incx,
                  &beta, y, incy);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintVectorZ(FULLPRINT, n, y, incy, "Y");

      mkl_free(a);
      mkl_free(x);
      mkl_free(y);

      return 0;
}

