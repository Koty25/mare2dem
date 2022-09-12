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
!      C B L A S _ S N R M 2   Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT  n, incx;
      float   *x;
      float    res;
      MKL_INT  len_x;

      printf("\n     C B L A S _ S N R M 2  EXAMPLE PROGRAM\n");

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
      if( GetIntegerParameters(in_file, &n, &incx) != 2 ) {
          printf("\n ERROR of n, incx reading\n");
          fclose(in_file);
          return 1;
      }

      len_x = 1+(n-1)*ABS(incx);
      x    = (float *)mkl_calloc(len_x, sizeof( float ), 64);
      if( x == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(x);
          fclose(in_file);
          return 1;
      }
      if( GetVectorS(in_file, n, x, 1) != len_x ) {
        printf("\n ERROR of vector X reading\n");
        mkl_free(x);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       N="INT_FORMAT, n);
      PrintVectorS(FULLPRINT, n, x, incx, "X");

/*      Call CBLAS_SNRM2 subroutine ( C Interface )            */

      res = cblas_snrm2(n, x, incx);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      printf("\n       SNRM2 = %7.3f", res);

      mkl_free(x);
      return 0;
}

