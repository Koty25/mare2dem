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
!      C B L A S _ S R O T M G   Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      float      param[5];
      float      sd1, sd2, sx1, sy1;

      param[0] = 0.0;
      param[1] = 0.0;
      param[2] = 0.0;
      param[3] = 0.0;
      param[4] = 0.0;

      printf("\n     C B L A S _ S R O T M G  EXAMPLE PROGRAM\n");

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
      if( GetScalarsS(in_file, &sd1, &sd2, &sx1, &sy1) != 4 ) {
          printf("\n ERROR of sd1, sd2, sx1, sy1 reading\n");
          fclose(in_file);
          return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       SD1=%5.2f  SD2=%5.2f  SX1=%5.2f  SY1=%5.2f",
               sd1, sd2, sx1, sy1);

/*      Call CBLAS_SROTMG subroutine ( C Interface )           */

      cblas_srotmg(&sd1, &sd2, &sx1, sy1, param);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      printf("\n       SD1=%5.2f  SD2=%5.2f  SX1=%5.2f",
               sd1, sd2, sx1);
      PrintVectorS(1, 5, param, 1, "PARAM");

      return 0;
}

