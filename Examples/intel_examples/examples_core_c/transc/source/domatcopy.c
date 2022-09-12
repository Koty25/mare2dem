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
!    mkl_domatcopy - out-of-place transposition routine,
!    Example Program Text ( C Interface )
!******************************************************************************/
#include <mkl_trans.h>
#include "common_func.h"

int main(int argc, char *argv[])
{ 
  size_t n=3, m=5;
  double src[] = { 
    1.,   2.,   3.,   4.,   5.,
    6.,   7.,   8.,   9.,   10.,
    11.,  12.,  13.,  14.,  15.
  }; /* source matrix */
  double dst[8]; /* destination matrix */
  size_t src_stride = 5;
  size_t dst_stride = 2;

  printf("\nThis is example of using mkl_domatcopy\n");

  printf("INPUT DATA:\nSource matrix:\n");
  print_matrix(n, m, 'd', src);

  /*
  **  Source submatrix(2,4) a will be transposed 
  */
  mkl_domatcopy('R'        /* row-major ordering */,
                'T'        /* A will be transposed */,
                2          /* rows */,
                4          /* cols */,
                1.         /* scales the input matrix */,
                src        /* source matrix */,
                src_stride /* src_stride */,
                dst        /* destination matrix */,
                dst_stride /* dst_stride */);
  /*  New matrix: src = { 
  **      1,  6,
  **      2,  7,
  **      3,  8,
  **      4,  9,
  **    }
  */
  printf("OUTPUT DATA:\nDestination matrix:\n");
  print_matrix(4, 2, 'd', dst);
  
  return 0;
}
