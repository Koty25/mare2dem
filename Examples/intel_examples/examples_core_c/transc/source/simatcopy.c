/*******************************************************************************
* Copyright 2010-2020 Intel Corporation.
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
!    mkl_simatcopy - in-place transposition routine,
!    Example Program Text ( C Interface )
!******************************************************************************/
#include <mkl_trans.h>
#include "common_func.h"

int main(int argc, char *argv[])
{
  size_t n=4, m=6; /* rows, cols of source matrix */
  float src[]= {                       
     1,  2,  3,  4,  5,  6,
     7,  8,  9,  10, 11, 12,
     13, 14, 15, 16, 17, 18,
     19, 20, 21, 22, 23, 24
  }; /* source matrix */

  printf("\nExample of using mkl_simatcopy transposition\n");
  printf("INPUT DATA:\n");
  printf("Source matrix:\n");
  print_matrix(n, m, 's', src);
                         
  /*
  **  Submatrix(3,4) will be transposed and the rest of the matrix will be unchanged    
  */
  mkl_simatcopy('R' /* row-major ordering */, 
                'T' /* A will be transposed */, 
                3   /* rows */, 
                4   /* cols */, 
                1.  /* scales the input matrix */, 
                src /* source matrix */, 
                6   /* src_lda */, 
                6   /* dst_lda */);
  /* New matrix: src = {  
  **    1,  7,  13, 4,  5,  6,
  **    2,  8,  14, 10, 11, 12,
  **    3,  9,  15, 16, 17, 18,
  **    4,  10, 16, 22, 23, 24        
  ** }
  */

  /*
  **  destination matrix will replace the source matrix with using dst_lda
  */
  printf("\nOUTPUT DATA:\nDestination matrix:\n");
  print_matrix(n, m, 's',src);

  return 0;

}
