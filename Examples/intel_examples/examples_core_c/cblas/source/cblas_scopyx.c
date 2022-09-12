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
!      C B L A S _ S C O P Y   Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT  n, incx, incy;
      float   *x, *y;
      MKL_INT  len_x, len_y;

      printf("\n     C B L A S _ S C O P Y  EXAMPLE PROGRAM\n");

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

      len_x = 1+(n-1)*ABS(incx);
      len_y = 1+(n-1)*ABS(incy);
      x    = (float *)mkl_calloc(len_x, sizeof( float ), 64);
      y    = (float *)mkl_calloc(len_y, sizeof( float ), 64);
      if( x == NULL || y == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(x);
          mkl_free(y);
          fclose(in_file);
          return 1;
      }
      if( GetVectorS(in_file, n, x, incx) != len_x ) {
        printf("\n ERROR of vector X reading\n");
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       N="INT_FORMAT, n);
      PrintVectorS(FULLPRINT, n, x, incx, "X");

/*      Call CBLAS_SCOPY subroutine ( C Interface )            */

      cblas_scopy(n, x, incx, y, incy);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintVectorS(FULLPRINT, n, y, incy, "Y");

      mkl_free(x);
      mkl_free(y);
      return 0;
}

