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
!    mkl_comatcopy2 - out-of-place transposition routine,
!    Example Program Text ( C Interface )
!******************************************************************************/
#include <mkl.h>
#include "common_func.h"

int main(int argc, char *argv[])
{
  printf("\nExample of using mkl_comatcopy2 transposition\n");

  size_t src_rows = 8; /* rows of source matrix */
  size_t src_cols = 8; /* cols of source matrix */
  size_t dst_rows = 4; /* rows of destination matrix */
  size_t dst_cols = 4; /* cols of destination matrix */
  MKL_Complex8 alpha;
  alpha.real = 1.;
  alpha.imag = 0.;
  size_t src_stride = 2; /* stride for source matrix */
  size_t dst_stride = 1; /* stride for destination matrix */

  /* Allocating and initializing source matrix */
  MKL_Complex8 *src = (MKL_Complex8*)mkl_malloc(src_rows*src_cols*sizeof(MKL_Complex8),64);
  for (int i = 0; i < src_rows; ++i){
    for (int j = 0; j < src_cols; ++j){
      src[i*src_cols + j].real = 10*(i+1) + (j+1);
      src[i*src_cols + j].imag = 10*(i+1) + (j+1);
    }
  }
  printf("INPUT DATA:\nSource matrix A:\n");
  print_matrix(src_rows, src_cols, 'c', src);
/*
    {11.,11.},   {12.,12.},   {13.,13.},   {14.,14.},   {15.,15.},   {16.,16.},   {17.,17.},   {18.,18.},
    {21.,21.},   {22.,22.},   {23.,23.},   {24.,24.},   {25.,25.},   {26.,26.},   {27.,27.},   {28.,28.},
    {31.,31.},   {32.,32.},   {33.,33.},   {34.,34.},   {35.,35.},   {36.,36.},   {37.,37.},   {38.,38.},
    {41.,41.},   {42.,42.},   {43.,43.},   {44.,44.},   {45.,45.},   {46.,46.},   {47.,47.},   {48.,48.},
    {51.,51.},   {52.,52.},   {53.,53.},   {54.,54.},   {55.,55.},   {56.,56.},   {57.,57.},   {58.,58.},
    {61.,61.},   {62.,62.},   {63.,63.},   {64.,64.},   {65.,65.},   {66.,66.},   {67.,67.},   {68.,68.},
    {71.,71.},   {72.,72.},   {73.,73.},   {74.,74.},   {75.,75.},   {76.,76.},   {77.,77.},   {78.,78.},
    {81.,81.},   {82.,82.},   {83.,83.},   {84.,84.},   {85.,85.},   {86.,86.},   {87.,87.},   {88.,88.}
*/

  MKL_Complex8 *dst = (MKL_Complex8*)mkl_malloc(dst_rows*dst_cols*sizeof(MKL_Complex8),64);

  printf("\nDestination matrix - elements from the source matrix are taken one by two each two rows with conjugate transposition:\n");
  mkl_comatcopy2('R'         /* row-major ordering */,
                 'C'         /* matrix A will be conjugate transposed */,
                  dst_rows   /* dst_rows */,
                  dst_cols   /* dst_cols */,
                  alpha      /* scales the input matrix by {1.,0.} */,
                  src        /* source matrix pointer */,
                  2*src_cols /* distance between adjacent rows in src */,
                  src_stride /* distance between adjacent columns in src */,
                  dst        /* destination matrix pointer */,
                  dst_cols   /* distance between adjacent rows in dst */,
                  dst_stride /* distance between adjacent columns in dst */);
  print_matrix(dst_cols, dst_rows, 'c', dst);
/*
    {11.,-11.},   {31.,-31.},   {51.,-51.},   {71.,-71.},
    {13.,-13.},   {33.,-33.},   {53.,-53.},   {73.,-73.},
    {15.,-15.},   {35.,-35.},   {55.,-55.},   {75.,-75.},
    {17.,-17.},   {37.,-37.},   {57.,-57.},   {77.,-77.}
*/
  printf("\nDestination matrix - elements from the source matrix are taken one by two each two rows with transposition (column-major):\n");
  mkl_comatcopy2('C'         /* column-major ordering */,
                 'T'         /* matrix A will be just transposed */,
                  dst_rows   /* dst_rows */,
                  dst_cols   /* dst_cols */,
                  alpha      /* scales the input matrix by {1.,0.} */,
                  src        /* source matrix pointer */,
                  src_stride /* distance between adjacent columns in src */,
                  2*src_cols /* distance between adjacent rows in src */,
                  dst        /* destination matrix pointer */,
                  dst_stride /* distance between adjacent columns in dst */,
                  dst_cols   /* distance between adjacent rows in dst */);
  print_matrix(dst_cols, dst_rows, 'c', dst);
  mkl_free(dst);
/*
    {11.,11.},   {31.,31.},   {51.,51.},   {71.,71.},
    {13.,13.},   {33.,33.},   {53.,53.},   {73.,73.},
    {15.,15.},   {35.,35.},   {55.,55.},   {75.,75.},
    {17.,17.},   {37.,37.},   {57.,57.},   {77.,77.}
*/

  dst_rows = 3;   /* rows of destination matrix */
  dst_cols = 10;  /* cols of destination matrix */
  alpha.real = -1.;
  alpha.imag = -1.;
  src_stride = 1; /* stride for source matrix */
  dst_stride = 2; /* stride for destination matrix */

  dst = (MKL_Complex8*)mkl_calloc((dst_rows+1)*dst_cols,sizeof(MKL_Complex8),64);

  printf("\nDestination matrix - copy of submatrix(3,5), first element is (4,4), putting with stride = 2 with conjugation:\n");
  mkl_comatcopy2('R'         /* row-major ordering */,
                 'R'         /* matrix A will be conjugated */,
                  dst_rows   /* dst_rows */,
                  dst_cols   /* dst_cols */,
                  alpha      /* scales the input matrix by {-1.,-1.} */,
                  src+3*src_cols+3 /* source matrix pointer, starting from element (4,4) */,
                  src_cols   /* distance between adjacent rows in src */,
                  src_stride /* distance between adjacent columns in src */,
                  dst        /* destination matrix pointer */,
                  dst_cols   /* distance between adjacent rows in dst */,
                  dst_stride /* distance between adjacent columns in dst */);
  print_matrix(dst_rows, dst_cols, 'c', dst);
  mkl_free(dst);
/*
    { -88.,0.},   {0.,0.},   { -90.,0.},   {0.,0.},   { -92.,0.},   {0.,0.},   { -94.,0.},   {0.,0.},
    {-108.,0.},   {0.,0.},   {-110.,0.},   {0.,0.},   {-112.,0.},   {0.,0.},   {-114.,0.},   {0.,0.},
    {-128.,0.},   {0.,0.},   {-130.,0.},   {0.,0.},   {-132.,0.},   {0.,0.},   {-134.,0.},   {0.,0.}
*/

#if defined(_OPENMP)
#include <omp.h>
  int num_of_threads = 2;
  dst_rows = 8;   /* rows of destination matrix */
  dst_cols = 8;   /* cols of destination matrix */
  alpha.real = 1.;
  alpha.imag = 0.;
  src_stride = 2; /* stride for source matrix */
  dst_stride = 2; /* stride for destination matrix */

  dst = (MKL_Complex8*)mkl_calloc((dst_rows+1)*dst_cols,sizeof(MKL_Complex8),64);
  printf("\nDestination matrix - same as original. Transposition parallel execution -\n");
  printf("    1st thread copies odd columns, 2nd one copies even columns:\n");
#pragma omp parallel num_threads(num_of_threads)
  {
    int thr_id = omp_get_thread_num();
    mkl_comatcopy2('R'         /* row-major ordering */,
                   'N'         /* matrix A will be simply copied */,
                    dst_rows   /* dst_rows */,
                    dst_cols   /* dst_cols */,
                    alpha      /* scales the input matrix by {1.,0.} */,
                    src+thr_id /* source matrix pointer, thread 0 starts from (0,0), thread 1 starts from (0,1) */,
                    src_cols   /* distance between adjacent rows in src */,
                    src_stride /* distance between adjacent columns in src */,
                    dst+thr_id /* destination matrix pointer, thread 0 starts (0,0), thread 1 starts from (0,1) */,
                    dst_cols   /* distance between adjacent rows in dst */,
                    dst_stride /* distance between adjacent columns in dst */);
  }
  print_matrix(dst_rows, dst_cols, 'c', dst);
  mkl_free(dst);
#endif

  mkl_free(src);
  return 0;
}
