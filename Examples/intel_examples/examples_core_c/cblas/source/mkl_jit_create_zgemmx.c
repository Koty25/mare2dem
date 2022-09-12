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
!      MKL_JIT_CREATE_ZGEMM
!      Example Program Text ( C Interface )
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#include "mkl_example.h"
#include "mkl_cblas.h"

int main(int argc, char *argv[])
{
      FILE *in_file;
      char *in_file_name;

      MKL_INT         m, n, k;
      MKL_INT         lda, ldb, ldc;
      MKL_INT         rmaxa, cmaxa, rmaxb, cmaxb, rmaxc, cmaxc;
      MKL_Complex16   alpha, beta;
      MKL_Complex16   *a, *b, *c;
      MKL_LAYOUT      layout;
      MKL_TRANSPOSE   transA, transB;
      MKL_INT         ma, na, mb, nb;
      void*           jitter;
      zgemm_jit_kernel_t  kernel;
      mkl_jit_status_t status;

      printf("\n     M K L _ J I T _ C R E A T E _ Z G E M M  EXAMPLE PROGRAM\n");

/*       Get input parameters                                  */

      if (argc == 1) {
          printf("\n Input data file must be specified as the first parameter");
          return 1;
      }
      in_file_name = argv[1];

/*       Get input data                                        */

      if ((in_file = fopen(in_file_name, "r")) == NULL) {
          printf("\n ERROR opening '%s' with mode=\"r\"\n", in_file_name);
          return 1;
      }
      if (GetIntegerParameters(in_file, &m, &n, &k) != 3) {
          printf("\n ERROR reading m, n, k\n");
          fclose(in_file);
          return 1;
      }
      if (GetScalarsZ(in_file, &alpha, &beta) != 2) {
          printf("\n ERROR reading alpha, beta\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, (CBLAS_TRANSPOSE*) &transA, (CBLAS_TRANSPOSE*) &transB, (CBLAS_LAYOUT*) &layout) != 3 ) {
          printf("\n ERROR reading transA, transB, layout\n");
          fclose(in_file);
          return 1;
      }

      if( transA == MKL_NOTRANS ) {
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
      if( transB == MKL_NOTRANS ) {
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
      a = (MKL_Complex16 *) mkl_calloc(rmaxa*cmaxa, sizeof(MKL_Complex16), 64);
      b = (MKL_Complex16 *) mkl_calloc(rmaxb*cmaxb, sizeof(MKL_Complex16), 64);
      c = (MKL_Complex16 *) mkl_calloc(rmaxc*cmaxc, sizeof(MKL_Complex16), 64);
      if(!a || !b || !c) {
          printf( "\n Can't allocate memory for arrays\n");
          mkl_free(a);
          mkl_free(b);
          mkl_free(c);
          fclose(in_file);
          return 1;
      }
      if( layout == MKL_ROW_MAJOR ) {
          lda=cmaxa;
          ldb=cmaxb;
          ldc=cmaxc;
      } else {
          lda=rmaxa;
          ldb=rmaxb;
          ldc=rmaxc;
      }
      if( GetArrayZ(in_file, (CBLAS_LAYOUT*) &layout, GENERAL_MATRIX, &ma, &na, a, &lda) != 0 ) {
        printf("\n ERROR reading array A\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayZ(in_file, (CBLAS_LAYOUT*) &layout, GENERAL_MATRIX, &mb, &nb, b, &ldb) != 0 ) {
        printf("\n ERROR reading array B\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      if( GetArrayZ(in_file, (CBLAS_LAYOUT*) &layout, GENERAL_MATRIX, &m,  &n,  c, &ldc) != 0 ) {
        printf("\n ERROR reading array C\n");
        mkl_free(a);
        mkl_free(b);
        mkl_free(c);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       M="INT_FORMAT"  N="INT_FORMAT"  K="INT_FORMAT, m, n, k);
      printf("\n       ALPHA=(%5.1f,%5.1f)  BETA=(%5.1f,%5.1f)", alpha.real, alpha.imag, beta.real, beta.imag);
      PrintParameters("TRANSA, TRANSB", transA, transB);
      PrintParameters("LAYOUT", layout);
      PrintArrayZ((CBLAS_LAYOUT*) &layout, FULLPRINT, GENERAL_MATRIX, &ma, &na, a, &lda, "A");
      PrintArrayZ((CBLAS_LAYOUT*) &layout, FULLPRINT, GENERAL_MATRIX, &mb, &nb, b, &ldb, "B");
      PrintArrayZ((CBLAS_LAYOUT*) &layout, FULLPRINT, GENERAL_MATRIX, &m,  &n,  c, &ldc, "C");

/*      Call MKL_JIT_CREATE_ZGEMM subroutine ( C Interface )                  */

      status = mkl_jit_create_zgemm(&jitter, layout, transA, transB, m, n, k, &alpha,
                                    lda, ldb, &beta, ldc);
      
      if (status == MKL_JIT_ERROR) {
          printf("Error: cannot create jitter\n");
          return 1;
      }

      kernel = mkl_jit_get_zgemm_ptr(jitter);

      kernel(jitter, a, b, c);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintArrayZ((CBLAS_LAYOUT*) &layout, FULLPRINT, GENERAL_MATRIX, &m, &n, c, &ldc, "C");

      mkl_jit_destroy(jitter);

      mkl_free(a);
      mkl_free(b);
      mkl_free(c);

      return 0;
}

