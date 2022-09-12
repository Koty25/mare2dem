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
!      C B L A S _ S A X P Y I  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT    n;
      float      alpha;
      float     *x, *y;
      MKL_INT   *indx;

      MKL_INT    indmax;

      printf("\n     C B L A S _ S A X P Y I  EXAMPLE PROGRAM\n");

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
      if( GetScalarsS(in_file, &alpha) != 1 ) {
          printf("\n ERROR of alpha reading\n");
          fclose(in_file);
          return 1;
      }

      x    = (float *)mkl_calloc(n, sizeof( float ), 64);
      indx = (MKL_INT *)mkl_calloc(n, sizeof( MKL_INT ), 64);
      if( x == NULL || indx == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(indx);
          mkl_free(x);
          fclose(in_file);
          return 1;
      }
      if( GetVectorS(in_file, n, x, 1) != n ) {
        printf("\n ERROR of vector X reading\n");
        mkl_free(indx);
        mkl_free(x);
        fclose(in_file);
        return 1;
      }
      if( GetVectorI(in_file, n, indx) != n ) {
        printf("\n ERROR of vector INDX reading\n");
        mkl_free(indx);
        mkl_free(x);
        fclose(in_file);
        return 1;
      }
      indmax = MaxValue(n, indx)+1;
      y    = (float *)mkl_calloc(indmax, sizeof( float ), 64);
      if( y == NULL ) {
          printf( "Can't allocate memory for arrays\n");
          mkl_free(indx);
          mkl_free(x);
          mkl_free(y);
          fclose(in_file);
          return 1;
      }
      if( GetVectorS(in_file, indmax, y, 1) != indmax ) {
        printf("\n ERROR of vector Y reading\n");
        mkl_free(indx);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       N="INT_FORMAT, n);
      printf("\n       ALPHA=%4.1f", alpha);
      PrintVectorS(SHORTPRINT, n, x, 1, "X");
      PrintVectorI(n, indx, "INDX");
      PrintVectorS(SHORTPRINT, indmax, y, 1, "Y");

/*      Call CBLAS_SAXPYI subroutine ( C Interface )           */

      cblas_saxpyi(n, alpha, x, indx, y);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintVectorS(SHORTPRINT, indmax, y, 1, "Y");

      mkl_free(x);
      mkl_free(indx);
      mkl_free(y);

      return 0;
}

