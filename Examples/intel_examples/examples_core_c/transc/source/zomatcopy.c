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
!    mkl_zomatcopy - out-of-place transposition routine,
!    Example Program Text ( C Interface )
!******************************************************************************/
#include <mkl_trans.h>
#include "common_func.h"

int main(int argc, char *argv[]) {  
  size_t n=4, m=6;
  MKL_Complex16 alpha;
  MKL_Complex16 src[] = {
      { 1.,  2.},  {3.,  4.},   {5.,  6.},  {7.,  8.},  {9., 10.}, {11., 12.},
      {13., 14.}, {15., 16.},  {17., 18.}, {19., 20.}, {21., 22.}, {23., 24.},
      {25., 26.}, {27., 28.},  {29., 30.}, {31., 32.}, {33., 34.}, {35., 36.},
      {37., 38.}, {39., 40.},  {41., 42.}, {43., 44.}, {45., 46.}, {47., 48.}
  };/* source matrix */
  MKL_Complex16 dst[12]; /* destination matrix */
  alpha.real = 1.;
  alpha.imag = 0.;

  printf("\nThis is example of using mkl_zomatcopy\n");
  printf("INPUT DATA:\nSource matrix:\n");
  print_matrix(n, m, 'z', src);

  /* 
  **  Source submatrix(3,4) a will be transposed
  */
  mkl_zomatcopy('R'   /* row-major ordering */, 
                'T'   /* A will be transposed */, 
                3     /* rows */, 
                4     /* cols */, 
                alpha /* scales the input matrix */, 
                src   /* source matrix */, 
                6     /* src_stride */, 
                dst   /* destination matrix */, 
                3     /* dst_stride */);
  /* New matrix: src = { 
  **    1., 2.,      13., 14.,    25., 26.,
  **    3., 4.,      15., 16.,    27., 28.,
  **    5., 6.,      17., 18.,    29., 30.,
  **    7., 8.,      19., 20.,    31., 32.
  **  }
  */
  printf("OUTPUT DATA:\nDestination matrix - source submatrix(3,4) a will be transposed:\n");
  print_matrix(4, 3, 'z', dst);

  /*
  **  Source submatrix(3,4) a will be conjugate transposed    
  */
  mkl_zomatcopy('R'   /* row-major ordering */, 
                'C'   /* A will be conjugate transposed */, 
                4     /* rows */, 
                3     /* cols */, 
                alpha /* scales the input matrix */, 
                src   /* source matrix */, 
                6     /* src_stride */, 
                dst   /* destination matrix */, 
                4     /* dst_stride */);
  /*  New matrix: src = { 
  **    1., -2.,      13., -14.,    25., -26.,  37., -38.,
  **    3., -4.,      15., -16.,    27., -28.,  39., -40.,
  **    5., -6.,      17., -18.,    29., -30.,  41., -42.
  **  }
  */
  printf("Destination matrix - source submatrix(4,3) a will be conjugate transposed:\n");
  print_matrix(3, 4, 'z', dst);

  /*
  **  Source submatrix(3,4) a will be conjugated    
  */
  mkl_zomatcopy('R'   /* row-major ordering */, 
                'R'   /* A will be conjugated */, 
                3     /* rows */, 
                4     /* cols */, 
                alpha /* scales the input matrix */, 
                src   /* source matrix */, 
                6     /* src_stride */, 
                dst   /* destination matrix */, 
                4     /* dst_stride */);
  /*  New matrix: src = {  
  **    1.,   -2.,   3.,  -4.,   5.,  -6.,   7.,  -8.,  
  **    13.,  -14., 15.,  -16., 17.,  -18., 19.,  -20., 
  **    25.,  -26., 27.,  -28., 29.,  -30., 31.,  -32.
  **  }
  */
  printf("Destination matrix - source submatrix(3,4) a will be conjugated:\n");
  print_matrix(3, 4, 'z', dst);
  
  return 0;
}
