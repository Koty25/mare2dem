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
!      C B L A S _ S R O T G   Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      float      c, s;
      float      sa, sb;

      printf("\n     C B L A S _ S R O T G  EXAMPLE PROGRAM\n");

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
      if( GetScalarsS(in_file, &sa, &sb) != 2 ) {
          printf("\n ERROR of sa, sb reading\n");
          fclose(in_file);
          return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       SA=%5.2f  SB=%5.2f", sa, sb);

/*      Call CBLAS_SROTG subroutine ( C Interface )            */

      cblas_srotg(&sa, &sb, &c, &s);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      printf("\n       SA=%6.3f  SB=%6.3f", sa, sb);
      printf("\n        C=%6.3f   S=%6.3f", c, s);

      return 0;
}

