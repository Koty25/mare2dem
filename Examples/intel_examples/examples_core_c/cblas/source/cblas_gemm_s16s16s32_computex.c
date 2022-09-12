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
!      C B L A S _ G E M M _ S 1 6 U 1 6 S 3 2 _ C O M P U T E  
!      Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT           m, n, k;
      MKL_INT           lda, ldb, ldc;
      MKL_INT           rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc, size_co;
      size_t            size_dest;
      float             alpha, beta;
      MKL_INT16         *a, *b, *dest;
      MKL_INT32         *c;
      CBLAS_LAYOUT      layout;
      CBLAS_IDENTIFIER  identifier;
      CBLAS_STORAGE     storage;
      CBLAS_TRANSPOSE   transA, transB;
      CBLAS_OFFSET      offsetc;
      MKL_INT           ma, na, mb, nb;
      MKL_INT16         ao, bo;
      MKL_INT32         *co;

      printf("\n     C B L A S _ G E M M _ S 1 6 U 1 6 S 3 2 _ C O M P U T E  EXAMPLE PROGRAM\n");

/*       Get input parameters                                                              */

      if( argc == 1 ) {
          printf("\n You must specify in_file data file as 1-st parameter");
          return 1;
      }
      in_file_name = argv[1];

/*       Get input data                                                                    */

      if( (in_file = fopen( in_file_name, "r" )) == NULL ) {
         printf("\n ERROR on OPEN '%s' with mode=\"r\"\n", in_file_name);
         return 1;
      }
      if( GetIntegerParameters(in_file, &m, &n, &k) != 3 ) {
          printf("\n ERROR of m, n, k reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetScalarsS(in_file, &alpha, &beta) != 2 ) {
          printf("\n ERROR of alpha, beta reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, &transA, &transB, &layout,  &offsetc) != 4 ) {
          printf("\n ERROR of transA, transB, layout, offsetc reading\n");
          fclose(in_file);
          return 1;
      }
      if (GetIntegerParameters16(in_file, &ao, &bo) != 2) {
          printf("\n ERROR of ao, bo reading\n");
          fclose(in_file);
          return 1;
      }

      if( transA == CblasNoTrans ) {
         rmaxa = m + 1;
         cmaxa = k;
         ma    = m;
         na    = k;
      } else {
         rmaxa = k + 1;
         cmaxa = m;
         ma    = k;
         na    = m;
      }
      if( transB == CblasNoTrans ) {
         rmaxb = k + 1;
         cmaxb = n;
         mb    = k;
         nb    = n;
      } else {
         rmaxb = n + 1;
         cmaxb = k;
         mb    = n;
         nb    = k;
      }
      rmaxc = m + 1;
      cmaxc = n;

      size_co = 1;
      if( offsetc == CblasRowOffset ) {
          size_co = n;
      }
      if( offsetc == CblasColOffset ) {
          size_co = m;
      }

      a = (MKL_INT16 *)mkl_calloc(rmaxa*cmaxa, sizeof( MKL_INT16 ), 64);
      b = (MKL_INT16 *)mkl_calloc(rmaxb*cmaxb, sizeof( MKL_INT16 ), 64);
      c = (MKL_INT32 *)mkl_calloc(rmaxc*cmaxc, sizeof( MKL_INT32 ), 64);
      co = (MKL_INT32 *)mkl_calloc(size_co, sizeof( MKL_INT32 ), 64);

      if( a == NULL || b == NULL || c == NULL || co == NULL ) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(a);
          mkl_free(b);
          mkl_free(c);
          mkl_free(co);
          fclose(in_file);
          return 1;
      }
      if (layout == CblasRowMajor) {
          lda=cmaxa;
          ldb=cmaxb;
          ldc=cmaxc;
      } else {
          lda=rmaxa;
          ldb=rmaxb;
          ldc=rmaxc;
      }
      if( GetVectorI32(in_file, size_co,  co) != size_co ) {
        printf("\n ERROR of vector co reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        mkl_free(co);
        fclose(in_file);
        return 1;
      }
      if( GetArrayI16(in_file, &layout, GENERAL_MATRIX, &ma, &na, a, &lda) != 0 ) {
        printf("\n ERROR of array A reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        mkl_free(co);
        fclose(in_file);
        return 1;
      }
      if( GetArrayI16(in_file, &layout, GENERAL_MATRIX, &mb, &nb, b, &ldb) != 0 ) {
        printf("\n ERROR of array B reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        mkl_free(co);
        fclose(in_file);
        return 1;
      }
      if( GetArrayI32(in_file, &layout, GENERAL_MATRIX, &m,  &n,  c, &ldc) != 0 ) {
        printf("\n ERROR of array C reading\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        mkl_free(co);
        fclose(in_file);
        return 1;
      }

      fclose(in_file);

/*       Print input data                                                                  */

      printf("\n     INPUT DATA");
      printf("\n       M="INT_FORMAT"  N="INT_FORMAT"  K="INT_FORMAT, m, n, k);
      printf("\n       ALPHA=%5.1f  BETA=%5.1f", alpha, beta);
      PrintParameters("TRANSA, TRANSB", transA, transB);
      PrintParameters("LAYOUT", layout);
      PrintParameters("OFFSETC", offsetc);
      printf("\n       ao=%d  bo=%d", (int) ao, (int) bo);
      PrintVectorI32(size_co, co, "co");
      PrintArrayI16(&layout, FULLPRINT, GENERAL_MATRIX, &ma, &na, a, &lda, "A");
      PrintArrayI16(&layout, FULLPRINT, GENERAL_MATRIX, &mb, &nb, b, &ldb, "B");
      PrintArrayI32(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, c, &ldc, "C");

/*       Call GEMM_S16U16S32_PACK_GET_SIZE function to allocate buffer for A (C Interface) */
      identifier = CblasAMatrix;
      size_dest = cblas_gemm_s16s16s32_pack_get_size(identifier, m, n, k);
      dest = (MKL_INT16 *)mkl_malloc(size_dest, 64);
      if( dest == NULL ) {
          printf( "\n Can't allocate memory for buffer\n");
          mkl_free(a);
          mkl_free(b);
          mkl_free(c);
          mkl_free(co);
          return 1;
      }

/*       Call GEMM_S16S16S32_PACK subroutine to perform packing (C Interface)              */
      cblas_gemm_s16s16s32_pack(layout, identifier, transA, m, n, k, a, lda, dest);

/*       Call GEMM_S16S16S32_COMPUTE subroutine (C Interface)                              */
      storage = CblasPacked;
      cblas_gemm_s16s16s32_compute(layout, storage, transB, offsetc, m, n, k,
                                   alpha, dest, lda, ao, b, ldb, bo, beta, c, ldc, co);

/*       Free allocated storage                                                            */
      mkl_free(dest);

/*       Print output data                                                                 */

      printf("\n\n     OUTPUT DATA");
      PrintArrayI32(&layout, FULLPRINT, GENERAL_MATRIX, &m, &n, c, &ldc, "C");

      mkl_free(a);
      mkl_free(b);
      mkl_free(c);
      mkl_free(co);

      return 0;
}

