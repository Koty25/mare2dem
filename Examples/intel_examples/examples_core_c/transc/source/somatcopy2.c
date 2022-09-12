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
!    mkl_somatcopy2 - out-of-place transposition routine,
!    Example Program Text ( C Interface )
!******************************************************************************/
#include <mkl.h>
#include "common_func.h"

int main(int argc, char *argv[])
{
  printf("\nExample of using mkl_domatcopy2 transposition\n");

  size_t src_rows = 8; /* rows of source matrix */
  size_t src_cols = 8; /* cols of source matrix */
  size_t dst_rows = 4; /* rows of destination matrix */
  size_t dst_cols = 4; /* cols of destination matrix */
  float alpha = 1.;
  size_t src_stride = 2; /* stride for source matrix */
  size_t dst_stride = 1; /* stride for destination matrix */

  /* Allocating and initializing source matrix */
  float *src = (float*)mkl_malloc(src_rows*src_cols*sizeof(float),64);
  for (int i = 0; i < src_rows; ++i){
    for (int j = 0; j < src_cols; ++j){
      src[i*src_cols + j] = 10*(i+1) + (j+1);
    }
  }
  printf("INPUT DATA:\nSource matrix A:\n");
  print_matrix(src_rows, src_cols, 's', src);
/*
    11,   12,   13,   14,   15,   16,   17,   18,
    21,   22,   23,   24,   25,   26,   27,   28,
    31,   32,   33,   34,   35,   36,   37,   38,
    41,   42,   43,   44,   45,   46,   47,   48,
    51,   52,   53,   54,   55,   56,   57,   58,
    61,   62,   63,   64,   65,   66,   67,   68,
    71,   72,   73,   74,   75,   76,   77,   78,
    81,   82,   83,   84,   85,   86,   87,   88,
*/

  float *dst = (float*)mkl_malloc(dst_rows*dst_cols*sizeof(float),64);

  printf("\nDestination matrix - elements from the source matrix are taken one by two each two rows with transposition:\n");
  mkl_somatcopy2('R'         /* row-major ordering */,
                 'T'         /* matrix A will be transposed */,
                  dst_rows   /* dst_rows */,
                  dst_cols   /* dst_cols */,
                  alpha      /* scales the input matrix by 1. */,
                  src        /* source matrix pointer */,
                  2*src_cols /* distance between adjacent rows in src */,
                  src_stride /* distance between adjacent columns in src */,
                  dst        /* destination matrix pointer */,
                  dst_cols   /* distance between adjacent rows in dst */,
                  dst_stride /* distance between adjacent columns in dst */);
  print_matrix(dst_cols, dst_rows, 's', dst);
/*
    11,   31,   51,   71,
    13,   33,   53,   73,
    15,   35,   55,   75,
    17,   37,   57,   77
*/
  printf("\nDestination matrix - same as previous but with column-major function call:\n");
  mkl_somatcopy2('C'         /* column-major ordering */,
                 'T'         /* matrix A will be transposed */,
                  dst_rows   /* dst_rows */,
                  dst_cols   /* dst_cols */,
                  alpha      /* scales the input matrix by 1. */,
                  src        /* source matrix pointer */,
                  src_stride /* distance between adjacent columns in src */,
                  2*src_cols /* distance between adjacent rows in src */,
                  dst        /* destination matrix pointer */,
                  dst_stride /* distance between adjacent columns in dst */,
                  dst_cols   /* distance between adjacent rows in dst */);
  print_matrix(dst_cols, dst_rows, 's', dst);
  mkl_free(dst);
/*
    11,   31,   51,   71,
    13,   33,   53,   73,
    15,   35,   55,   75,
    17,   37,   57,   77
*/

  dst_rows = 3;   /* rows of destination matrix */
  dst_cols = 10;  /* cols of destination matrix */
  alpha = -1.;
  src_stride = 1; /* stride for source matrix */
  dst_stride = 2; /* stride for destination matrix */

  dst = (float*)mkl_calloc((dst_rows+1)*dst_cols,sizeof(float),64);

  printf("\nDestination matrix - copy of submatrix(3,5), first element is (4,4), putting with stride = 2:\n");
#if 1
  mkl_somatcopy2('R'         /* row-major ordering */,
                 'N'         /* matrix A will be simply copied */,
                  dst_rows   /* dst_rows */,
                  dst_cols   /* dst_cols */,
                  alpha      /* scales the input matrix by -1. */,
                  src+3*src_cols+3 /* source matrix pointer, starting from element (4,4) */,
                  src_cols   /* distance between adjacent rows in src */,
                  src_stride /* distance between adjacent columns in src */,
                  dst        /* destination matrix pointer */,
                  dst_cols   /* distance between adjacent rows in dst */,
                  dst_stride /* distance between adjacent columns in dst */);
#endif
  print_matrix(dst_rows, dst_cols, 's', dst);
  mkl_free(dst);
/*
    -44,   0,   -45,   0,   -46,   0,   -47,   0,   -48,   0,
    -54,   0,   -55,   0,   -56,   0,   -57,   0,   -58,   0,
    -64,   0,   -65,   0,   -66,   0,   -67,   0,   -68,   0
*/


#if defined(_OPENMP)
#include <omp.h>
  int num_of_threads = 2;
  dst_rows = 8;   /* rows of destination matrix */
  dst_cols = 8;   /* cols of destination matrix */
  alpha = 1.;
  src_stride = 2; /* stride for source matrix */
  dst_stride = 2; /* stride for destination matrix */

  dst = (float*)mkl_calloc((dst_rows+1)*dst_cols,sizeof(float),64);
  printf("\nDestination matrix - same as original. Transposition parallel execution -\n");
  printf("    1st thread copies odd columns, 2nd one copies even columns:\n");
#pragma omp parallel num_threads(num_of_threads)
  {
    int thr_id = omp_get_thread_num();
    mkl_somatcopy2('R'         /* row-major ordering */,
                   'N'         /* matrix A will be simply copied */,
                    dst_rows   /* dst_rows */,
                    dst_cols   /* dst_cols */,
                    alpha      /* scales the input matrix by 1. */,
                    src+thr_id /* source matrix pointer, thread 0 starts from (0,0), thread 1 starts from (0,1) */,
                    src_cols   /* distance between adjacent rows in src */,
                    src_stride /* distance between adjacent columns in src */,
                    dst+thr_id /* destination matrix pointer, thread 0 starts (0,0), thread 1 starts from (0,1) */,
                    dst_cols   /* distance between adjacent rows in dst */,
                    dst_stride /* distance between adjacent columns in dst */);
  }
  print_matrix(dst_rows, dst_cols, 's', dst);
//  mkl_free(dst);
#endif
  mkl_free(src);
  return 0;
}
