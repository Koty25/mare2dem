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
! Content:
!   mkl_cimatcopy - in-place transposition routine,
!   Example Program Text ( C Interface )
!******************************************************************************/
#include <mkl_trans.h>
#include "common_func.h"

int main(int argc, char *argv[])
{
  size_t n=4, m=6; /* rows, cols of source matrix */
  MKL_Complex8 alpha;
  MKL_Complex8 src[] = {
      { 1.,  2.},  {3.,  4.},   {5.,  6.},  {7.,  8.},  {9., 10.}, {11., 12.},
      {13., 14.}, {15., 16.},  {17., 18.}, {19., 20.}, {21., 22.}, {23., 24.},
      {25., 26.}, {27., 28.},  {29., 30.}, {31., 32.}, {33., 34.}, {35., 36.},
      {37., 38.}, {39., 40.},  {41., 42.}, {43., 44.}, {45., 46.}, {47., 48.}
  };/* source matrix */
  alpha.real = 1.;
  alpha.imag = 0.;

  printf("\nExample of using mkl_cimatcopy transposition\n");
  printf("INPUT DATA:\n");
  printf("Source matrix:\n");
  print_matrix(n, m, 'c', src);
                       
  /* 
  **  Submatrix(3,4) will be transposed and the rest of the matrix will be unchanged    
  */
  mkl_cimatcopy('R'    /* row-major ordering */, 
                'T'    /* matrix will be transposed */, 
                 3     /* rows */, 
                 4     /* cols */, 
                 alpha /* scales the input matrix */, 
                 src   /* source matrix */, 
                 6     /* src_lda */, 
                 6     /* dst_lda */);
  /* New matrix: src = {  
  **   1.,  2., 13., 14., 25., 26.,  7.,  8.,  9., 10., 11., 12.,
  **   3.,  4., 15., 16., 27., 28., 19., 20., 21., 22., 23., 24.,
  **   5.,  6., 17., 18., 29., 30., 31., 32., 33., 34., 35., 36.,
  **   7.,  8., 19., 20., 31., 32., 43., 44., 45., 46., 47., 48.      
  ** }
  */
  printf("OUTPUT DATA:\nDestination matrix - submatrix(3,4) will be transposed and the rest of the matrix will be unchanged:\n");
  print_matrix(n, m, 'c', src);

  /*
  **  Submatrix(2,6) will be conjugated and the rest of the matrix will be unchanged    
  */
  mkl_cimatcopy('R'    /* row-major ordering */, 
                'R'    /* matrix will be conjugated */, 
                 2     /* rows */, 
                 6     /* cols */, 
                 alpha /* scales the input matrix */, 
                 src   /* source matrix */, 
                 6     /* src_lda */, 
                 6     /* dst_lda */);
  /* New matrix: src = { 
  **    1., -2.,  13., -14.,  25., -26.,   7.,  -8.,   9., -10.,  11., -12.,
  **    3., -4.,  15., -16.,  27., -28.,  19.,  -20., 21., -22.,  23., -24.,
  **    5.,  6.,  17.,  18.,  29.,  30.,  31.,   32., 33.,  34.,  35.,  36.,
  **    7.,  8.,  19.,  20.,  31.,  32.,  43.,   44., 45.,  46.,  47.,  48.     
  **  }
  */
  /*  
  **  Destination matrix will replace the source matrix with using dst_lda
  */
  printf("Destination matrix - submatrix(2,6) will be conjugated and the rest of the matrix will be unchanged:\n");
  print_matrix(n, m, 'c', src);

  /*
  **  Submatrix(2,6) will be conjugate transposed and the rest of the matrix will be unchanged    
  */
  mkl_cimatcopy('R'    /* row-major ordering */, 
                'C'    /* matrix will be conjugate transposed */, 
                 2     /* rows */, 
                 3     /* cols */, 
                 alpha /* scales the input matrix */, 
                 src   /* source matrix */, 
                 6     /* src_lda */, 
                 6     /* dst_lda */);
  /* New matrix: src = { 
  **    1.,  2.,   3.,  4., 25., -26.,   7.,  -8.,   9., -10.,  11., -12.,
  **   13.,  14., 15., 16., 27., -28.,  19., -20.,  21., -22.,  23., -24.,
  **   25.,  26., 27., 28., 29.,  30.,  31.,  32.,  33.,  34.,  35.,  36.,
  **    7.,  8.,  19., 20., 31.,  32.,  43.,  44.,  45.,  46.,  47.,  48.     
  **  }
  */
  /* 
  **  Destination matrix will replace the source matrix with using dst_lda
  */
  printf("Destination matrix - submatrix(2,6) will be conjugate transposed and the rest of the matrix will be unchanged:\n");
  print_matrix(n, m, 'c', src);

  return 0;
}
