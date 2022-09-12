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
!      C B L A S _ D R O T M G   Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      double     param[5];
      double     dd1, dd2, dx1, dy1;

      param[0] = 0.0;
      param[1] = 0.0;
      param[2] = 0.0;
      param[3] = 0.0;
      param[4] = 0.0;

      printf("\n     C B L A S _ D R O T M G  EXAMPLE PROGRAM\n");

/*       Get input parameters                                  */

      if( argc == 1 ) {
         printf("\n You must specify in_file data file as 1-st parameter");
         return 1;
      }
      in_file_name = argv[1];

/*       Get input data                                       */

      if( (in_file = fopen( in_file_name, "r" )) == NULL ) {
         printf("\n ERROR on OPEN '%s' with mode=\"r\"\n", in_file_name);
         return 1;
      }
      if( GetScalarsD(in_file, &dd1, &dd2, &dx1, &dy1) != 4 ) {
          printf("\n ERROR of dd1, dd2, dx1, dy1 reading\n");
          fclose(in_file);
          return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       DD1=%5.2f  DD2=%5.2f  DX1=%5.2f  DY1=%5.2f",
               dd1, dd2, dx1, dy1);

/*      Call CBLAS_DROTMG subroutine ( C Interface )            */

      cblas_drotmg(&dd1, &dd2, &dx1, dy1, param);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      printf("\n       DD1=%6.3f  DD2=%6.3f  DX1=%6.3f",  dd1, dd2, dx1);
      PrintVectorD(SHORTPRINT, 5, param, 1, "PARAM");

      return 0;
}

