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
!  mkl_domatadd - out-of-place transposition routine,
!  Example Program Text ( C Interface )
!******************************************************************************/
#include <mkl_trans.h>
#include "common_func.h"

int main(int argc, char *argv[])
{ 
  size_t n=3, m=5;/* rows, cols of source matrix  */
  double a[] = {
    1.,  2.,   3.,    4.,    50., 
    5.,  6.,   7.,    8.,    60.,
    9.,  10.,  11.,  12.,    70.
  }; /* source matrix A */
  double b[] = {     
    1.,  2.,   3.,    4.,    50., 
    5.,  6.,   7.,    8.,    40.,
    9.,  10.,  11.,  12.,    30.
  }; /* source matrix B */
  double dst[9]; /* destination matrix */

  printf("\nExample of using mkl_domatadd transposition\n");
  printf("INPUT DATA:\nSource matrix A:\n");
  print_matrix(n, m, 'd', a);

  printf("Source matrix B:\n");
  print_matrix(n, m, 'd', b);

  /*
  **  Addition of transposed sub-matrix(3,3) a and unchanged sub-matrix(3,3) b
  */
  mkl_domatadd('R'  /* row-major ordering */, 
               'T'  /* A will be transposed */, 
               'N'  /* no changes to B */, 
                3   /* rows */, 
                3   /* cols */, 
                1.  /* alpha */, 
                a   /* source matrix */, 
                5   /* lda */, 
                1.  /* beta */, 
                b   /* source matrix */, 
                5   /* ldb */, 
                dst /* destination matrix */, 
                3   /* ldc */); 
  /* New matrix: c =  { 
  **    2.,  7.,   12.,
  **    7.,  12.,  17.,
  **    12., 17.,  22.
  ** }
  */
  printf("OUTPUT DATA:\nDestination matrix - addition of transposed submatrix(3,3) of A and submatrix of B:\n");  
  print_matrix(3, 3, 'd', dst);
  /*
  **  Addition of transposed sub-matrices(3,3) a and b
  */
  mkl_domatadd('R'  /* row-major ordering */, 
               'T'  /* A will be transposed */, 
               'T'  /* B will be transposed */, 
                3   /* rows */, 
                3   /* cols */, 
                1.  /* alpha */, 
                a   /* source matrix */, 
                5   /* lda */, 
                1.  /* beta */, 
                b   /* source matrix */, 
                5   /* ldb */, 
                dst /* destination matrix */, 
                3   /* ldc */); 
  /* New matrix: c =  { 
  **    2., 10.,  18.,
  **    4., 12.,  20.,
  **    6., 14.,  22.
  **    }
  */
  printf("Destination matrix - Addition of transposed submatrices(3,3) of A and B:\n"); 
  print_matrix(3, 3, 'd', dst);
  
  return 0;
}
