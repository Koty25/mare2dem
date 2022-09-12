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
!      C B L A S _ C D O T C I  Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT       n;
      MKL_INT      *indx;
      MKL_Complex8 *x, *y;

      MKL_INT       indmax;
      MKL_Complex8  res;

      res.real = 0.;
      res.imag = 0.;

      printf("\n     C B L A S _ C D O T C I  EXAMPLE PROGRAM\n");

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

      x    = (MKL_Complex8 *)mkl_calloc(n, sizeof( MKL_Complex8 ), 64);
      indx = (MKL_INT *)mkl_calloc(n, sizeof( MKL_INT ), 64);
      if( x == NULL || indx == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(indx);
          mkl_free(x);
          fclose(in_file);
          return 1;
      }
      if( GetVectorC(in_file, n, x, 1) != n ) {
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
      y    = (MKL_Complex8 *)mkl_calloc(indmax, sizeof( MKL_Complex8 ), 64);
      if( y == NULL ) {
          printf( "Can't allocate memory for arrays\n");
          mkl_free(indx);
          mkl_free(x);
          mkl_free(y);
          fclose(in_file);
          return 1;
      }
      if( GetVectorC(in_file, indmax, y, 1) != indmax ) {
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
      PrintVectorC(SHORTPRINT, n, x, 1, "X");
      PrintVectorI(n, indx, "INDX");
      PrintVectorC(SHORTPRINT, indmax, y, 1, "Y");

/*      Call CBLAS_CDOTCI subroutine ( C Interface )           */

      cblas_cdotci_sub(n, x, indx, y, &res);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      printf("\n       CDOTCI = (%6.2f,%6.2f)", res.real, res.imag);

      mkl_free(x);
      mkl_free(indx);
      mkl_free(y);

      return 0;
}

