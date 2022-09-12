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
!  mkl_comatadd - out-of-place transposition routine,
!  Example Program Text ( C Interface )
!******************************************************************************/
#include <mkl_trans.h>
#include "common_func.h"

int main(int argc, char *argv[])
{ 
  size_t n=4, m=3; /* rows, cols of source matrix */
  MKL_Complex8 alpha;
  MKL_Complex8 beta;
  MKL_Complex8 a[]={
    { 1.,  2.},  { 3.,  4.}, { 5.,  6.},
    { 7.,  8.},  { 9., 10.}, {11., 12.},
    {13., 14.},  {15., 16.}, {17., 18.},
    {25., 26.},  {27., 28.}, {29., 30.}
  }; /* source matrix */
  MKL_Complex8 b[]={
    { 1.1,  2.1}, { 3.2,  4.2}, { 5.3,  6.3},
    { 7.1,  8.1}, { 9.2, 10.2}, {11.3, 12.3},
    {13.1, 14.1}, {15.2, 16.2}, {17.3, 18.3},
    {25.1, 26.1}, {27.2, 28.2}, {29.3, 30.3}
  }; /* source matrix   */
  MKL_Complex8 dst[9];/* destination matrix */
  alpha.real = 1.;
  alpha.imag = 0.;
  beta.real = 1.;
  beta.imag = 0.;

  printf("\nExample of using mkl_comatadd transposition\n");
  printf("INPUT DATA:\nSource matrix A:\n");
  print_matrix(n, m, 'c', a);

  printf("Source matrix B:\n");
  print_matrix(n, m, 'c', b);

  /* 
  **  Addition of transposed submatrix(3,3) a and unchanged submatrix(3,3) b 
  */
  mkl_comatadd('R'    /* row-major ordering */, 
               'T'    /* A will be transposed */,
               'N'    /* no changes to B */, 
                3     /* rows */, 
                3     /* cols */, 
                alpha /* alpha */, 
                a     /* source matrix */, 
                3     /* lda */, 
                beta  /* beta */, 
                b     /* source matrix */, 
                3     /* ldb */, 
                dst   /* destination matrix */, 
                3     /* ldc */); 
  /*
  **  New matrix: c =  { 
  **    2.1, 4.1,     10.2, 12.2,    18.3, 20.3,
  **    10.1, 12.1,    18.2, 20.2,    26.3, 28.3,
  **    18.1, 20.1,    26.2, 28.2,    34.3, 36.3
  **  }
  */
  printf("OUTPUT DATA:\nDestination matrix - addition of transposed submatrix(3,3) A and submatrix B:\n");  
  print_matrix(3, 3, 'c', dst);
  /*
  **  Addition of transposed submatrix(3,3) a and conjugate transposed submatrix(3,3) b 
  */
  mkl_comatadd('R'    /* row-major ordering */, 
               'T'    /* A will be transposed */,
               'C'    /* B will be conjugate transposed */, 
                3     /* rows */, 
                3     /* cols */, 
                alpha /* alpha */, 
                a     /* source matrix */, 
                3     /* lda */, 
                beta  /* beta */, 
                b     /* source matrix */, 
                3     /* ldb */, 
                dst   /* destination matrix */, 
                3     /* ldc */);  
  /*  New matrix: c = {
  **    2.1, -0.1,  14.1, -0.1,  26.1, -0.1,
  **    6.2, -0.2,  18.2, -0.2,  30.2, -0.2,
  **    10.3, -0.3, 22.3, -0.3,  34.3, -0.3   
  **  }
  */
  printf("Destination matrix - Addition of transposed submatrix(3,3) A and conjugate transposed submatrix(3,3) B:\n");  
  print_matrix(3, 3, 'c', dst);
  
  return 0;
}
