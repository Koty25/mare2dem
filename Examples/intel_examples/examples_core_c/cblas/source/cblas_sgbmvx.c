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
!      C B L A S _ S G B M V  Example Program Text ( C Interface )
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

      MKL_INT         m, n, kl, ku, lda, incx, incy;
      MKL_INT         rmaxa, cmaxa;
      float           alpha, beta;
      float          *a, *x, *y;
      CBLAS_LAYOUT    layout;
      CBLAS_TRANSPOSE trans;
      MKL_INT         nx, ny, len_x, len_y;

      printf("\n     C B L A S _ S G B M V  EXAMPLE PROGRAM\n");

/*       Get input parameters                                  */

      if( argc == 1 ) {
         printf("\n You must specify in_file data file as 1-st parameter");
         return 1;
      }
      in_file_name = argv[1];

/*       Get input data                                        */

      if( (in_file = fopen( in_file_name, "r" )) == NULL ) {
         printf("\n ERROR on OPEN '%s' with mode=\"r\"\n", in_file_name);
         return 1;
      }
      if( GetIntegerParameters(in_file, &m, &n, &kl, &ku) != 4 ) {
          printf("\n ERROR of m, n, kl, ku reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetIntegerParameters(in_file, &incx, &incy) != 2 ) {
          printf("\n ERROR of incx, incy reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetScalarsS(in_file, &alpha, &beta) != 2 ) {
          printf("\n ERROR of alpha, beta reading\n");
          fclose(in_file);
          return 1;
      }
      if( GetCblasCharParameters(in_file, &trans, &layout) != 2 ) {
          printf("\n ERROR of trans, layout reading\n");
          fclose(in_file);
          return 1;
      }

      rmaxa = kl + ku + 1;
      if( layout == CblasRowMajor )
          cmaxa = m;
      else
          cmaxa = n;
      a = (float *)mkl_calloc(rmaxa*cmaxa, sizeof(float), 64);
      if( trans == CblasNoTrans ) {
         nx = n;
         ny = m;
      } else {
         nx = m;
         ny = n;
      }
      len_x = 1+(nx-1)*ABS(incx);
      len_y = 1+(ny-1)*ABS(incy);
      x = (float *)mkl_calloc(len_x, sizeof(float), 64);
      y = (float *)mkl_calloc(len_y, sizeof(float), 64);
      if( a == NULL || x == NULL || y == NULL ) {
          printf("\n Can't allocate memory for arrays\n");
          mkl_free(a);
          mkl_free(x);
          mkl_free(y);
          fclose(in_file);
          return 1;
      }
      if( GetVectorS(in_file, nx, x, incx) != len_x ) {
        printf("\n ERROR of vector X reading \n");
        mkl_free(a);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      if( GetVectorS(in_file, ny, y, incy) != len_y ) {
        printf("\n ERROR of vector Y reading\n");
        mkl_free(a);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      lda=rmaxa;
      if( GetBandArrayS(in_file, &layout, kl, ku, m, n, a, lda) != 0 ) {
        printf("\n ERROR of array A reading\n");
        mkl_free(a);
        mkl_free(x);
        mkl_free(y);
        fclose(in_file);
        return 1;
      }
      fclose(in_file);

/*       Print input data                                      */

      printf("\n     INPUT DATA");
      printf("\n       M="INT_FORMAT"  N="INT_FORMAT, m, n);
      printf("\n       ALPHA =%4.1f  BETA =%4.1f", alpha, beta);
      PrintParameters("TRANS", trans);
      PrintParameters("LAYOUT", layout);
      PrintVectorS(FULLPRINT, nx, x, incx, "X");
      PrintVectorS(FULLPRINT, ny, y, incy, "Y");
      PrintBandArrayS(&layout, FULLPRINT, kl, ku, m, n, a, lda, "A");

/*      Call CBLAS_SGBMV subroutine ( C Interface )            */

      cblas_sgbmv(layout, trans, m, n, kl, ku, alpha, a, lda,
                  x, incx, beta, y, incy);

/*       Print output data                                     */

      printf("\n\n     OUTPUT DATA");
      PrintVectorS(FULLPRINT, ny, y, incy, "Y");

      mkl_free(a);
      mkl_free(x);
      mkl_free(y);

      return 0;
}

