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
!      C B L A S _ D S P M V  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT         n, incx, incy;
      MKL_INT         apmax;
      double          alpha, beta;
      double         *ap, *x, *y;
      CBLAS_LAYOUT    layout;
      CBLAS_UPLO      uplo;
      MKL_INT         len_x, len_y;

      printf("\n     C B L A S _ D S P M V  EXAMPLE PROGRAM\n");

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
      if( GetIntegerParameters(in_file, &incx, &incy) != 2 ) {
          printf("\n ERROR of incx, incy reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetScalarsD(in_file, &alpha, &beta) != 2 ) {
          printf("\n ERROR of alpha, beta reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, &uplo, &layout) != 2 ) {
          printf("\n ERROR of uplo, layout reading\n");
          fclose(in_file);
          return 1;
      }

      apmax = (n*(n+1))/2;
      len_x = 1+(n-1)*ABS(incx);
      len_y = 1+(n-1)*ABS(incy);
      ap = (double *)mkl_calloc(apmax, sizeof(double), 64);
      x  = (double *)mkl_calloc(len_x, sizeof(double), 64);
      y  = (double *)mkl_calloc(len_y, sizeof(double), 64);
      if( ap == NULL || x == NULL || y == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(ap);
          mkl_free(x);
          mkl_free(y);
          fclose(in_file);
          return 1;
      }
      if( GetVectorD(in_file, n, x, incx) != len_x ) {
        printf("\n ERROR of vector X reading\n");
        mkl_free(ap);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      if( GetVectorD(in_file, n, y, incy) != len_y ) {
        printf("\n ERROR of vector Y reading\n");
        mkl_free(ap);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      if( GetVectorD(in_file, (n*(n+1))/2, ap, 1) != (n*(n+1))/2 ) {
        printf("\n ERROR of vector AP reading\n");
        mkl_free(ap);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       N="INT_FORMAT, n);
      printf("\n       ALPHA=%5.1f  BETA=%5.1f", alpha, beta);
      PrintParameters("UPLO", uplo);
      PrintParameters("LAYOUT", layout);
      PrintVectorD(FULLPRINT, n, x, incx, "X");
      PrintVectorD(FULLPRINT, n, y, incy, "Y");
      PrintVectorD(SHORTPRINT, (n*(n+1))/2, ap, 1, "AP");

/*      Call CBLAS_DSPMV subroutine ( C Interface )            */

      cblas_dspmv(layout, uplo, n, alpha, ap, x, incx, beta, y, incy);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintVectorD(FULLPRINT, n, y, incy, "Y");

      mkl_free(ap);
      mkl_free(x);
      mkl_free(y);

      return 0;
}

