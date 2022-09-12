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
!
!******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <math.h>

#include "mkl_example.h"
#include "mkl_cblas.h"

MKL_INT MaxValue(MKL_INT n, MKL_INT *x)
{
      MKL_INT   i, indmax;

      indmax = x[0];
      for (i = 1; i < n; i++)
          if (indmax < x[i]) indmax = x[i];
      return indmax;
} /* MaxValue */

typedef union {
  float float_part;
  MKL_BF16 int_part[2];
} conv_t;

static float b2f(MKL_BF16 src) {
  conv_t conv;
  conv.int_part[0] = 0;
  conv.int_part[1] = src;
  return conv.float_part;
}

static MKL_BF16 f2b(float src) {
  conv_t conv;
  conv.float_part = src;
  return conv.int_part[1];
}

MKL_INT GetVectorI(FILE *in_file, MKL_INT n, MKL_INT *x)
{
      return GetValuesI( in_file, x, 0, n );
} /* GetVectorI */

MKL_INT GetVectorI16(FILE *in_file, MKL_INT n, MKL_INT16 *x)
{
      return GetValuesI16( in_file, x, 1, 0, n );
} /* GetVectorI16 */

MKL_INT GetVectorI32(FILE *in_file, MKL_INT n, MKL_INT32 *x)
{
      return GetValuesI32( in_file, x, 1, 0, n );
} /* GetVectorI32 */

MKL_INT GetVectorS(FILE *in_file, MKL_INT n, float *x,  MKL_INT incx)
{
      return GetValuesS( in_file, x, 1, 0, (1+(n-1)*ABS(incx)) );
} /* GetVectorS */

MKL_INT GetVectorD(FILE *in_file, MKL_INT n, double *x,  MKL_INT incx)
{
      return GetValuesD( in_file, x, 1, 0, (1+(n-1)*ABS(incx)) );
} /* GetVectorD */

MKL_INT GetVectorC(FILE *in_file, MKL_INT n, MKL_Complex8 *x,  MKL_INT incx)
{
      return GetValuesC( in_file, x, 1, 0, (1+(n-1)*ABS(incx)) );
} /* GetVectorC */

MKL_INT GetVectorZ(FILE *in_file, MKL_INT n, MKL_Complex16 *x,  MKL_INT incx)
{
      return GetValuesZ( in_file, x, 1, 0, (1+(n-1)*ABS(incx)) );
} /* GetVectorZ */


MKL_INT GetArrayI8(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT flag, MKL_INT *m, MKL_INT *n,
                    MKL_INT8 *a, MKL_INT *lda)
{
      MKL_INT i, j , number;
      MKL_INT8 *addr;

      if( *layout == CblasRowMajor ) {
        if( flag == GENERAL_MATRIX ) {
           for( i = 0; i < (*m); i++ ) {
              addr = a + i*(*lda);
              number = GetValuesI8( in_file, addr, 1, 0, *n );
              if( number != *n ) return 1;
           } /* for */
        } else if( flag == UPPER_MATRIX ) {
           for (i = 0; i < (*m); i++) {
              addr = a + i*(*lda);
              number = GetValuesI8( in_file, addr, 1, i, *n-i );
              if( number != *n-i ) return 1;
           } /* for */
        } else if( flag == LOWER_MATRIX ) {
           for (i = 0; i < (*m); i++) {
              addr = a + i*(*lda);
              number = GetValuesI8( in_file, addr, 1, 0, i+1 );
              if( number != i+1 ) return 1;
           } /* for */
        } /* if */
     } else if( *layout == CblasColMajor ) {
        if( flag == GENERAL_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI8( in_file, addr, 1, 0, *m );
              if( number != *m ) return 1;
           } /* for */
        } else if( flag == UPPER_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI8( in_file, addr, 1, 0, j+1 );
              if( number != j+1 ) return 1;
           } /* for */
        } else if( flag == LOWER_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI8( in_file, addr, 1, j, *m-j );
              if( number != *m-j ) return 1;
           } /* for */
        } /* if */
     } /* if */
     return 0;
} /* GetArrayI8 */

MKL_INT GetArrayI16(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT flag, MKL_INT *m, MKL_INT *n,
                    MKL_INT16 *a, MKL_INT *lda)
{
      MKL_INT i, j , number;
      MKL_INT16 *addr;

      if( *layout == CblasRowMajor ) {
        if( flag == GENERAL_MATRIX ) {
           for( i = 0; i < (*m); i++ ) {
              addr = a + i*(*lda);
              number = GetValuesI16( in_file, addr, 1, 0, *n );
              if( number != *n ) return 1;
           } /* for */
        } else if( flag == UPPER_MATRIX ) {
           for (i = 0; i < (*m); i++) {
              addr = a + i*(*lda);
              number = GetValuesI16( in_file, addr, 1, i, *n-i );
              if( number != *n-i ) return 1;
           } /* for */
        } else if( flag == LOWER_MATRIX ) {
           for (i = 0; i < (*m); i++) {
              addr = a + i*(*lda);
              number = GetValuesI16( in_file, addr, 1, 0, i+1 );
              if( number != i+1 ) return 1;
           } /* for */
        } /* if */
     } else if( *layout == CblasColMajor ) {
        if( flag == GENERAL_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI16( in_file, addr, 1, 0, *m );
              if( number != *m ) return 1;
           } /* for */
        } else if( flag == UPPER_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI16( in_file, addr, 1, 0, j+1 );
              if( number != j+1 ) return 1;
           } /* for */
        } else if( flag == LOWER_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI16( in_file, addr, 1, j, *m-j );
              if( number != *m-j ) return 1;
           } /* for */
        } /* if */
     } /* if */
     return 0;
} /* GetArrayI16 */

MKL_INT GetArrayI32(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT flag, MKL_INT *m, MKL_INT *n,
                    MKL_INT32 *a, MKL_INT *lda)
{
      MKL_INT i, j , number;
      MKL_INT32 *addr;

      if( *layout == CblasRowMajor ) {
        if( flag == GENERAL_MATRIX ) {
           for( i = 0; i < (*m); i++ ) {
              addr = a + i*(*lda);
              number = GetValuesI32( in_file, addr, 1, 0, *n );
              if( number != *n ) return 1;
           } /* for */
        } else if( flag == UPPER_MATRIX ) {
           for (i = 0; i < (*m); i++) {
              addr = a + i*(*lda);
              number = GetValuesI32( in_file, addr, 1, i, *n-i );
              if( number != *n-i ) return 1;
           } /* for */
        } else if( flag == LOWER_MATRIX ) {
           for (i = 0; i < (*m); i++) {
              addr = a + i*(*lda);
              number = GetValuesI32( in_file, addr, 1, 0, i+1 );
              if( number != i+1 ) return 1;
           } /* for */
        } /* if */
     } else if( *layout == CblasColMajor ) {
        if( flag == GENERAL_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI32( in_file, addr, 1, 0, *m );
              if( number != *m ) return 1;
           } /* for */
        } else if( flag == UPPER_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI32( in_file, addr, 1, 0, j+1 );
              if( number != j+1 ) return 1;
           } /* for */
        } else if( flag == LOWER_MATRIX ) {
           for( j = 0; j < (*n); j++ ) {
              addr = a + j*(*lda);
              number = GetValuesI32( in_file, addr, 1, j, *m-j );
              if( number != *m-j ) return 1;
           } /* for */
        } /* if */
     } /* if */
     return 0;
} /* GetArrayI32 */

MKL_INT GetArrayBF16(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT flag, MKL_INT *m, MKL_INT *n,
              MKL_BF16 *a,  MKL_INT *lda)
{
      MKL_INT   i, j, number;
      MKL_BF16 *addr;

      if( *layout == CblasRowMajor ) {
         if( flag == GENERAL_MATRIX ) {
            for( i = 0; i < (*m); i++ ) {
               addr = a + i*(*lda);
               number = GetValuesBF16( in_file, addr, 1, 0, *n );
               if( number != *n ) return 1;
            } /* for */
         } else if( flag == UPPER_MATRIX ) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesBF16( in_file, addr, 1, i, *n-i );
               if( number != *n-i ) return 1;
            } /* for */
         } else if( flag == LOWER_MATRIX ) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesBF16( in_file, addr, 1, 0, i+1 );
               if( number != i+1 ) return 1;
            } /* for */
         } /* if */
      } else if( *layout == CblasColMajor ) {
         if( flag == GENERAL_MATRIX ) {
            for( j = 0; j < (*n); j++ ) {
               addr = a + j*(*lda);
               number = GetValuesBF16( in_file, addr, 1, 0, *m );
               if( number != *m ) return 1;
            } /* for */
         } else if( flag == UPPER_MATRIX ) {
            for( j = 0; j < (*n); j++ ) {
               addr = a + j*(*lda);
               number = GetValuesBF16( in_file, addr, 1, 0, j+1 );
               if( number != j+1 ) return 1;
            } /* for */
         } else if( flag == LOWER_MATRIX ) {
            for( j = 0; j < (*n); j++ ) {
               addr = a + j*(*lda);
               number = GetValuesBF16( in_file, addr, 1, j, *m-j );
               if( number != *m-j ) return 1;
            } /* for */
         } /* if */
      } /* if */
      return 0;
} /* GetArrayBF16 */

MKL_INT GetArrayS(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT flag, MKL_INT *m, MKL_INT *n,
              float *a,  MKL_INT *lda)
{
      MKL_INT   i, j, number;
      float    *addr;

      if( *layout == CblasRowMajor ) {
         if( flag == GENERAL_MATRIX ) {
            for( i = 0; i < (*m); i++ ) {
               addr = a + i*(*lda);
               number = GetValuesS( in_file, addr, 1, 0, *n );
               if( number != *n ) return 1;
            } /* for */
         } else if( flag == UPPER_MATRIX ) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesS( in_file, addr, 1, i, *n-i );
               if( number != *n-i ) return 1;
            } /* for */
         } else if( flag == LOWER_MATRIX ) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesS( in_file, addr, 1, 0, i+1 );
               if( number != i+1 ) return 1;
            } /* for */
         } /* if */
      } else if( *layout == CblasColMajor ) {
         if( flag == GENERAL_MATRIX ) {
            for( j = 0; j < (*n); j++ ) {
               addr = a + j*(*lda);
               number = GetValuesS( in_file, addr, 1, 0, *m );
               if( number != *m ) return 1;
            } /* for */
         } else if( flag == UPPER_MATRIX ) {
            for( j = 0; j < (*n); j++ ) {
               addr = a + j*(*lda);
               number = GetValuesS( in_file, addr, 1, 0, j+1 );
               if( number != j+1 ) return 1;
            } /* for */
         } else if( flag == LOWER_MATRIX ) {
            for( j = 0; j < (*n); j++ ) {
               addr = a + j*(*lda);
               number = GetValuesS( in_file, addr, 1, j, *m-j );
               if( number != *m-j ) return 1;
            } /* for */
         } /* if */
      } /* if */
      return 0;
} /* GetArrayS */

MKL_INT GetArrayD(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT flag, MKL_INT *m, MKL_INT *n,
              double *a,  MKL_INT *lda)
{
      MKL_INT   i, j, number;
      double   *addr;

      if (*layout == CblasRowMajor) {
         if (flag == GENERAL_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesD( in_file, addr, 1, 0, *n );
               if( number != *n ) return 1;
            } /* for */
         } else if (flag == UPPER_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesD( in_file, addr, 1, i, *n-i );
               if( number != *n-i ) return 1;
            } /* for */
         } else if (flag == LOWER_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesD( in_file, addr, 1, 0, i+1 );
               if( number != i+1 ) return 1;
            } /* for */
         } /* if */
      } else if (*layout == CblasColMajor) {
         if (flag == GENERAL_MATRIX) {
             for (j = 0; j < (*n); j++) {
                addr = a + j*(*lda);
                number = GetValuesD( in_file, addr, 1, 0, *m );
                if( number != *m ) return 1;
             } /* for */
         } else if (flag == UPPER_MATRIX) {
             for (j = 0; j < (*n); j++) {
                 addr = a + j*(*lda);
                 number = GetValuesD( in_file, addr, 1, 0, j+1 );
                 if( number != j+1 ) return 1;
             } /* for */
         } else if (flag == LOWER_MATRIX) {
             for (j = 0; j < (*n); j++) {
                 addr = a + j*(*lda);
                 number = GetValuesD( in_file, addr, 1, j, *m-j );
                 if( number != *m-j ) return 1;
             } /* for */
         } /* if */
      } /* if */
      return 0;
} /* GetArrayD */

MKL_INT GetArrayC(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT flag, MKL_INT *m, MKL_INT *n,
              MKL_Complex8 *a,  MKL_INT *lda)
{
      MKL_INT        i, j, number;
      MKL_Complex8  *addr;

      if (*layout == CblasRowMajor) {
         if (flag == GENERAL_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesC( in_file, addr, 1, 0, *n );
               if( number != *n ) return 1;
            } /* for */
         } else if (flag == UPPER_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesC( in_file, addr, 1, i, *n-i );
               if( number != *n-i ) return 1;
            } /* for */
         } else if (flag == LOWER_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesC( in_file, addr, 1, 0, i+1 );
               if( number != i+1 ) return 1;
            } /* for */
         } /* if */
      } else if (*layout == CblasColMajor) {
         if (flag == GENERAL_MATRIX) {
            for (j = 0; j < (*n); j++) {
               addr = a + j*(*lda);
               number = GetValuesC( in_file, addr, 1, 0, *m );
               if( number != *m ) return 1;
            } /* for */
         } else if (flag == UPPER_MATRIX) {
            for (j = 0; j < (*n); j++) {
               addr = a + j*(*lda);
               number = GetValuesC( in_file, addr, 1, 0, j+1 );
               if( number != j+1 ) return 1;
            } /* for */
         } else if (flag == LOWER_MATRIX) {
            for (j = 0; j < (*n); j++) {
               addr = a + j*(*lda);
               number = GetValuesC( in_file, addr, 1, j, *m-j );
               if( number != *m-j ) return 1;
            } /* for */
         } /* if */
      } /* if */
      return 0;
} /* GetArrayC */

MKL_INT GetArrayZ(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT flag,
              MKL_INT *m, MKL_INT *n, MKL_Complex16 *a,  MKL_INT *lda)
{
      MKL_INT        i, j, number;
      MKL_Complex16 *addr;

      if (*layout == CblasRowMajor) {
         if (flag == GENERAL_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesZ( in_file, addr, 1, 0, *n );
               if( number != *n ) return 1;
            } /* for */
         } else if (flag == UPPER_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesZ( in_file, addr, 1, i, *n-i );
               if( number != *n-i ) return 1;
            } /* for */
         } else if (flag == LOWER_MATRIX) {
            for (i = 0; i < (*m); i++) {
               addr = a + i*(*lda);
               number = GetValuesZ( in_file, addr, 1, 0, i+1 );
               if( number != i+1 ) return 1;
            } /* for */
         } /* if */
      } else if(*layout == CblasColMajor) {
         if (flag == GENERAL_MATRIX) {
            for (j = 0; j < (*n); j++) {
               addr = a + j*(*lda);
               number = GetValuesZ( in_file, addr, 1, 0, *m );
               if( number != *m ) return 1;
            } /* for */
         } else if (flag == UPPER_MATRIX) {
            for (j = 0; j < (*n); j++) {
               addr = a + j*(*lda);
               number = GetValuesZ( in_file, addr, 1, 0, j+1 );
               if( number != j+1 ) return 1;
            } /* for */
         } else if (flag == LOWER_MATRIX) {
            for (j = 0; j < (*n); j++) {
               addr = a + j*(*lda);
               number = GetValuesZ( in_file, addr, 1, j, *m-j );
               if( number != *m-j ) return 1;
            } /* for */
         } /* if */
      } /* if */
      return 0;
} /* GetArrayZ */

MKL_INT GetBandArrayS(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT kl, MKL_INT ku,
                  MKL_INT m, MKL_INT n, float *a,  MKL_INT lda)
{
      MKL_INT  i, j;
      MKL_INT  kl1, ku1, i_start, j_start, j_end, ku_rows, kl_rows, number;
      float   *addr, *addr1;

      if (*layout == CblasRowMajor) {
         for( i = 0; i < MIN( m, n ); i++ ) {
            addr = a + i*lda;
            kl1  = ( i - kl <= 0 ) ? i : kl;
            ku1  = ( i + ku >= n ) ? MAX(0,n-i-1) : ku;
            j_start = kl - kl1;
            j_end = j_start + kl1 + ku1;
            addr1 = addr + j_start;
            number = GetValuesS( in_file, addr1, 1, 0, j_end-j_start+1 );
            if( number != j_end-j_start+1 ) return 1;
         }
         for( i = MIN( m, n ); i < MIN( m, MIN( m, n ) + kl); i++ ) {
            addr = a + i*lda;
            kl1  = n - i + kl;
            j_start = ( kl > n ) ? kl - n : n - kl;
            j_end = j_start + kl1 - 1;
            addr1 = addr + j_start;
            number = GetValuesS( in_file, addr1, 1, 0, j_end-j_start+1 );
            if( number != j_end-j_start+1 ) return 1;
         }
      } else if (*layout == CblasColMajor) {
         i_start = (ku > n ) ? ku - n + 1 : 0;
         ku_rows = (ku > n ) ? n - 1 : ku;
         j_start = ku_rows;
         for( i = 0; i < ku_rows; i++ ) {
              j = j_start*lda; addr1 = a+i+i_start;
              number = GetValuesS( in_file, addr1, lda, j, n-ku_rows );
              if( number != n-ku_rows ) return 1;
              j_start--;
         }

         j_end = MIN(m,n);
         addr1 = a + ku;
         number = GetValuesS( in_file, addr1, lda, 0, j_end );
         if( number != j_end ) return 1;

         kl_rows = ( kl <= m-1 ) ?  kl : m - 1;
         for ( i = 1; i < kl_rows+1; i++ ) {
              kl1 = ( i+j_end <= m ) ? j_end : m - i;
              addr1 = a + ku + i;
              number = GetValuesS( in_file, addr1, lda, 0, kl1 );
              if( number != kl1 ) return 1;
         }
      }
      return 0;
} /* GetBandArrayS */

MKL_INT GetBandArrayD(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT kl, MKL_INT ku,
                  MKL_INT m, MKL_INT n, double *a,  MKL_INT lda)
{
      MKL_INT     i, j;
      MKL_INT     kl1, ku1, i_start, j_start, j_end, kl_rows, ku_rows, number;
      double     *addr, *addr1;

      if (*layout == CblasRowMajor) {
         for( i = 0; i < MIN( m, n ); i++ ) {
            addr = a + i*lda;
            kl1  = ( i - kl <= 0 ) ? i : kl;
            ku1  = ( i + ku >= n ) ? MAX(0,n-i-1) : ku;
            j_start = kl - kl1;
            j_end = j_start + kl1 + ku1;
            addr1 = addr + j_start;
            number = GetValuesD( in_file, addr1, 1, 0, j_end-j_start+1 );
            if( number != j_end-j_start+1 ) return 1;
         }
         for( i = MIN( m, n ); i < MIN( m, MIN( m, n ) + kl); i++ ) {
            addr = a + i*lda;
            kl1  = n - i + kl;
            j_start = ( kl > n ) ? kl - n : n - kl;
            j_end = j_start + kl1 - 1;
            addr1 = addr + j_start;
            number = GetValuesD( in_file, addr1, 1, 0, j_end-j_start+1 );
            if( number != j_end-j_start+1 ) return 1;
         }
      } else if (*layout == CblasColMajor) {
         i_start = (ku > n ) ? ku - n + 1 : 0;
         ku_rows = (ku > n ) ? n - 1 : ku;
         j_start = ku_rows;
         for( i = 0; i < ku_rows; i++ ) {
              j = j_start*lda; addr1 = a+i+i_start;
              number = GetValuesD( in_file, addr1, lda, j, n-ku_rows );
              if( number != n-ku_rows ) return 1;
              j_start--;
         }

         j_end = MIN(m,n);
         addr1 = a+ku;
         number = GetValuesD( in_file, addr1, lda, 0, j_end );
         if( number != j_end ) return 1;

         kl_rows = ( kl <= m-1 ) ?  kl : m - 1;
         for ( i = 1; i < kl_rows+1; i++ ) {
              kl1 = ( i+j_end <= m ) ? j_end : m - i;
              addr1 = a+ku+i;
              number = GetValuesD( in_file, addr1, lda, 0, kl1 );
              if( number != kl1 ) return 1;
         }
      } /* if */
      return 0;
} /* GetBandArrayD */

MKL_INT GetBandArrayC(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT kl, MKL_INT ku,
                  MKL_INT m, MKL_INT n, MKL_Complex8 *a,  MKL_INT lda)
{
      MKL_INT        i, j;
      MKL_INT        kl1, ku1, i_start, j_start, j_end, ku_rows, kl_rows, number;
      MKL_Complex8  *addr, *addr1;

      if (*layout == CblasRowMajor) {
         for( i = 0; i < MIN( m, n ); i++ ) {
            addr = a + i*lda;
            kl1  = ( i - kl <= 0 ) ? i : kl;
            ku1  = ( i + ku >= n ) ? MAX(0,n-i-1) : ku;
            j_start = kl - kl1;
            j_end = j_start + kl1 + ku1;
            addr1 = addr + j_start;
            number = GetValuesC( in_file, addr1, 1, 0, j_end-j_start+1 );
            if( number != j_end-j_start+1 ) return 1;
         }
         for( i = MIN( m, n ); i < MIN( m, MIN( m, n ) + kl); i++ ) {
            addr = a + i*lda;
            kl1  = n - i + kl;
            j_start = ( kl > n ) ? kl - n : n - kl;
            j_end = j_start + kl1 - 1;
            addr1 = addr + j_start;
            number = GetValuesC( in_file, addr1, 1, 0, j_end-j_start+1 );
            if( number != j_end-j_start+1 ) return 1;
         }
      } else if (*layout == CblasColMajor) {
         i_start = (ku > n ) ? ku - n + 1 : 0;
         ku_rows = (ku > n ) ? n - 1 : ku;
         j_start = ku_rows;
         for( i = 0; i < ku_rows; i++ ) {
              j = j_start*lda; addr1 = a+i+i_start;
              number = GetValuesC( in_file, addr1, lda, j, n-ku_rows );
              if( number != n-ku_rows ) return 1;
              j_start--;
         }

         j_end = MIN(m,n);
         addr1 = a+ku;
         number = GetValuesC( in_file, addr1, lda, 0, j_end );
         if( number != j_end ) return 1;

         kl_rows = ( kl <= m-1 ) ?  kl : m - 1;
         for ( i = 1; i < kl_rows+1; i++ ) {
              kl1 = ( i+j_end <= m ) ? j_end : m - i;
              addr1 = a+ku+i;
              number = GetValuesC( in_file, addr1, lda, 0, kl1 );
              if( number != kl1 ) return 1;
         }
      } /* if */
      return 0;
} /* GetBandArrayC */

MKL_INT GetBandArrayZ(FILE *in_file, CBLAS_LAYOUT*layout, MKL_INT kl, MKL_INT ku,
                  MKL_INT m, MKL_INT n, MKL_Complex16 *a,  MKL_INT lda)
{
      MKL_INT        i, j;
      MKL_INT        kl1, ku1, i_start, j_start, j_end, ku_rows, kl_rows, number;
      MKL_Complex16 *addr, *addr1;

      if (*layout == CblasRowMajor) {
         for( i = 0; i < MIN( m, n ); i++ ) {
            addr = a + i*lda;
            kl1  = ( i - kl <= 0 ) ? i : kl;
            ku1  = ( i + ku >= n ) ? MAX(0,n-i-1) : ku;
            j_start = kl - kl1;
            j_end = j_start + kl1 + ku1;
            addr1 = addr + j_start;
            number = GetValuesZ( in_file, addr1, 1, 0, j_end-j_start+1 );
            if( number != j_end-j_start+1 ) return 1;
         }
         for( i = MIN( m, n ); i < MIN( m, MIN( m, n ) + kl); i++ ) {
            addr = a + i*lda;
            kl1  = n - i + kl;
            j_start = ( kl > n ) ? kl - n : n - kl;
            j_end = j_start + kl1 - 1;
            addr1 = addr + j_start;
            number = GetValuesZ( in_file, addr1, 1, 0, j_end-j_start+1 );
            if( number != j_end-j_start+1 ) return 1;
         }
      } else if (*layout == CblasColMajor) {
         i_start = (ku > n ) ? ku - n + 1 : 0;
         ku_rows = (ku > n ) ? n - 1 : ku;
         j_start = ku_rows;
         for( i = 0; i < ku_rows; i++ ) {
              j = j_start*lda; addr1 = a+i+i_start;
              number = GetValuesZ( in_file, addr1, lda, j, n-ku_rows );
              if( number != n-ku_rows ) return 1;
              j_start--;
         }

         j_end = MIN(m,n);
         addr1 = a+ku;
         number = GetValuesZ( in_file, addr1, lda, 0, j_end );
         if( number != j_end ) return 1;

         kl_rows = ( kl <= m-1 ) ?  kl : m - 1;
         for ( i = 1; i < kl_rows+1; i++ ) {
              kl1 = ( i+j_end <= m ) ? j_end : m - i;
              addr1 = a+ku+i;
              number = GetValuesZ( in_file, addr1, lda, 0, kl1 );
              if( number != kl1 ) return 1;
         }
      } /* if */
      return 0;
} /* GetBandArrayZ */

MKL_INT GetValuesI( FILE *in_file, MKL_INT *in_array, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    int     value;
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%d", &value) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i]=(MKL_INT)value;
       counter++;
       if ( (str = strtok( NULL, " " )) == NULL  ) break;
    }
    return counter;
}

MKL_INT GetValuesI8( FILE *in_file, MKL_INT8 *in_array, MKL_INT ld, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    int value_int;
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%d", &value_int) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i*ld]=(MKL_INT8)value_int;
       counter++;
       if ( (str = strtok( NULL, " " )) == NULL  ) break;
    }
    return counter;
}



MKL_INT GetValuesI16( FILE *in_file, MKL_INT16 *in_array, MKL_INT ld, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    MKL_INT16     value;
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%hd", &value) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i*ld]=(MKL_INT16)value;
       counter++;
       if ( (str = strtok( NULL, " " )) == NULL  ) break;
    }
    return counter;
}

MKL_INT GetValuesI32( FILE *in_file, MKL_INT32 *in_array, MKL_INT ld, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    MKL_INT32     value;
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%d", &value) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i*ld]=(MKL_INT32)value;
       counter++;
       if ( (str = strtok( NULL, " " )) == NULL  ) break;
    }
    return counter;
}

MKL_INT GetValuesBF16( FILE *in_file, MKL_BF16 *in_array, MKL_INT ld, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    float   temp;
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%f", &temp) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i*ld]=f2b(temp);
       counter++;
       if ( (str = strtok( NULL, " " )) == NULL  ) break;
    }
    return counter;
} /* GetValuesBF16 */

MKL_INT GetValuesD( FILE *in_file, double *in_array, MKL_INT ld, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    double  temp;
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%lf", &temp) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i*ld]=temp;
       counter++;
       if ( (str = strtok( NULL, " " )) == NULL  ) break;
    }
    return counter;
} /* GetValuesD */

MKL_INT GetValuesS( FILE *in_file, float *in_array, MKL_INT ld, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    float   temp;
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%f", &temp) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i*ld]=temp;
       counter++;
       if ( (str = strtok( NULL, " " )) == NULL  ) break;
    }
    return counter;
} /* GetValuesS */

MKL_INT GetValuesC( FILE *in_file, MKL_Complex8 *in_array, MKL_INT ld, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    float   re, im;
    char    seps[]=" ,()\t";
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, seps );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%f", &re) != 1 || (str = strtok( NULL,seps)) == NULL ||
           sscanf( str, "%f", &im) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i*ld].real=re;
       in_array[begin+i*ld].imag=im;
       counter++;
       if ( (str = strtok( NULL, seps )) == NULL  ) break;
    }
    return counter;
} /* GetValuesC */

MKL_INT GetValuesZ( FILE *in_file, MKL_Complex16 *in_array, MKL_INT ld, MKL_INT begin, MKL_INT max_numbers )
{
    MKL_INT i, counter=0;
    double  re, im;
    char    seps[]=" ,()\t";
    char    buf[MAX_STRING_LEN], *str;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, seps );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    for( i = 0; i < max_numbers; i++ ) {
       if ( *str==COMMENTS ) break;
       if( sscanf( str, "%lf", &re) != 1 || (str = strtok( NULL,seps))== NULL ||
           sscanf( str, "%lf", &im) != 1 ){
           printf( "\n File format is inappropriate\n" );
           return 0;
       }
       in_array[begin+i*ld].real=re;
       in_array[begin+i*ld].imag=im;
       counter++;
       if ( (str = strtok( NULL, seps )) == NULL  ) break;
    }
    return counter;
} /* GetValuesZ */

int GetIntegerParameters8(FILE *in_file, ...)
{
    int       counter=0;
    char      buf[MAX_STRING_LEN], *str;
    MKL_INT8 *param;
    va_list   ap;
    int      value;

    do {
       fgets(buf, MAX_STRING_LEN, in_file);
       str = strtok( buf, " " );
       if( str == NULL ){
           printf("\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start(ap, in_file);
    param = va_arg(ap, MKL_INT8 *);

	while ( 1 ) {
        if( *str == COMMENTS ) break;
        if( sscanf(str, "%d", &value ) != 1 ){
            printf("\n File format is inappropriate\n");
            return 0;
        }
        *param = (MKL_INT8)value;
        counter++;
        if( (str = strtok( NULL, " " )) == NULL ) break;
        param = va_arg( ap, MKL_INT8 * );
    }
    va_end( ap );
    return counter;
} /* GetIntegerParameters16 */

int GetIntegerParameters16(FILE *in_file, ...)
{
    int       counter=0;
    char      buf[MAX_STRING_LEN], *str;
    MKL_INT16 *param;
    va_list   ap;
    MKL_INT16      value;

    do {
       fgets(buf, MAX_STRING_LEN, in_file);
       str = strtok( buf, " " );
       if( str == NULL ){
           printf("\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start(ap, in_file);
    param = va_arg(ap, MKL_INT16 *);

	while ( 1 ) {
        if( *str == COMMENTS ) break;
        if( sscanf(str, "%hd", &value ) != 1 ){
            printf("\n File format is inappropriate\n");
            return 0;
        }
        *param = (MKL_INT16)value;
        counter++;
        if( (str = strtok( NULL, " " )) == NULL ) break;
        param = va_arg( ap, MKL_INT16 * );
    }
    va_end( ap );
    return counter;
} /* GetIntegerParameters16 */

int GetIntegerParameters32(FILE *in_file, ...)
{
    int       counter=0;
    char      buf[MAX_STRING_LEN], *str;
    MKL_INT32 *param;
    va_list   ap;
    MKL_INT32      value;

    do {
       fgets(buf, MAX_STRING_LEN, in_file);
       str = strtok( buf, " " );
       if( str == NULL ){
           printf("\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start(ap, in_file);
    param = va_arg(ap, MKL_INT32 *);

	while ( 1 ) {
        if( *str == COMMENTS ) break;
        if( sscanf(str, "%d", &value ) != 1 ){
            printf("\n File format is inappropriate\n");
            return 0;
        }
        *param = (MKL_INT32)value;
        counter++;
        if( (str = strtok( NULL, " " )) == NULL ) break;
        param = va_arg( ap, MKL_INT32 * );
    }
    va_end( ap );
    return counter;
} /* GetIntegerParameters32 */

int GetIntegerParameters(FILE *in_file, ...)
{
    int      counter=0;
    char     buf[MAX_STRING_LEN], *str;
    MKL_INT *param;
    va_list  ap;
    int      value;

    do {
       fgets(buf, MAX_STRING_LEN, in_file);
       str = strtok( buf, " " );
       if( str == NULL ){
           printf("\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start(ap, in_file);
    param = va_arg(ap, MKL_INT *);

	while ( 1 ) {
        if( *str == COMMENTS ) break;
        if( sscanf(str, "%d", &value ) != 1 ){
            printf("\n File format is inappropriate\n");
            return 0;
        }
        *param = (MKL_INT)value;
        counter++;
        if( (str = strtok( NULL, " " )) == NULL ) break;
        param = va_arg( ap, MKL_INT * );
    }
    va_end( ap );
    return counter;
} /* GetIntegerParameters */

int GetCblasCharParameters(FILE *in_file, ...)
{
    int      counter=0;
    char     buf[MAX_STRING_LEN], *str;
    int     *param;
    va_list  ap;

    do {
       fgets(buf, MAX_STRING_LEN, in_file);
       str = strtok( buf, " " );
       if( str == NULL ){
           printf("\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start(ap, in_file);
    param = va_arg(ap, int *);

	while ( 1 ) {
        if( *str == COMMENTS ) break;
        if( sscanf(str, "%d", param ) != 1 ){
            printf("\n File format is inappropriate\n");
            return 0;
        }
        counter++;
        if( (str = strtok( NULL, " " )) == NULL ) break;
        param = va_arg( ap, int * );
    }
    va_end( ap );
    return counter;
} /* GetCblasCharParameters */

int GetScalarsD( FILE *in_file, ... )
{
    int     counter=0;
    char    buf[MAX_STRING_LEN], *str;
    double  *param;
    va_list ap;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start( ap, in_file );
    param = va_arg( ap, double * );
    while ( 1 ) {
        if ( *str == COMMENTS ) break;
        if( sscanf( str, "%lf", param ) != 1 ){
            printf( "\n File format is inappropriate\n" );
            return 0;
        }
        counter++;
        if ( (str = strtok( NULL, " " )) == NULL  ) break;
        param = va_arg( ap, double * );
    }
    va_end( ap );
    return counter;
} /* GetScalarsD */

int GetScalarsS( FILE *in_file, ... )
{
    int     counter=0;
    char    buf[MAX_STRING_LEN], *str;
    float  *param;
    va_list ap;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, " " );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start( ap, in_file );
    param = va_arg( ap, float * );
    while ( 1 ) {
        if ( *str == COMMENTS )  break;
        if( sscanf( str, "%f", param ) != 1 ){
            printf( "\n File format is inappropriate\n" );
            return 0;
        }
        counter++;
        if ( (str = strtok( NULL, " " )) == NULL  ) break;
        param = va_arg( ap, float * );
    }
    va_end( ap );
    return counter;
} /* GetScalarsS */

int GetScalarsC( FILE *in_file, ... )
{
    int           counter=0;
    char          buf[MAX_STRING_LEN], *str;
    char          seps[]=" ,()\t";
    float         re, im;
    MKL_Complex8 *param;
    va_list       ap;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, seps );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start( ap, in_file );
    param = va_arg( ap, MKL_Complex8 * );
    while ( 1 ) {
        if ( *str == COMMENTS )  break;
        if( sscanf( str, "%f", &re) != 1 || (str = strtok( NULL,seps))== NULL ||
            sscanf( str, "%f", &im) != 1 ){
            printf( "\n File format is inappropriate %f %f\n", re, im );
            return 0;
        }
        param->real = re;
        param->imag = im;
        counter++;
        if ( (str = strtok( NULL, seps )) == NULL  )  break;
        param = va_arg( ap, MKL_Complex8 * );
    }
    va_end( ap );
    return counter;
} /* GetScalarsC */

int GetScalarsZ( FILE *in_file, ... )
{
    int            counter=0;
    char           buf[MAX_STRING_LEN], *str;
    char           seps[]=" ,()\t";
    double         re, im;
    MKL_Complex16 *param;
    va_list        ap;

    do {
       fgets( buf, MAX_STRING_LEN, in_file );
       str = strtok( buf, seps );
       if( str == NULL ){
           printf( "\n File format is inappropriate\n");
           return 0;
       }
    } while ( *str == COMMENTS );
    va_start( ap, in_file );
    param = va_arg( ap, MKL_Complex16 * );
    while ( 1 ) {
        if ( *str == COMMENTS ) break;
        if( sscanf( str, "%lf", &re) != 1 || (str = strtok( NULL,seps)) == NULL ||
            sscanf( str, "%lf", &im) != 1 ){
            printf( "\n File format is inappropriate %f %f\n", re, im );
            return 0;
        }
        param->real = re;
        param->imag = im;
        counter++;
        if ( (str = strtok( NULL, seps )) == NULL  ) break;
        param = va_arg( ap, MKL_Complex16 * );
    }
    va_end( ap );
    return counter;
} /* GetScalarsZ */

void PrintParameters( char *names, ... )
{
    char            *p, *str, str1[MAX_STRING_LEN], buf[MAX_STRING_LEN];
    va_list          ap;
    char             seps[]=" ,";
    MKL_INT          itmp;
    MKL_UINT         i;

    printf("\n       ");
    va_start( ap, names );
    strcpy(buf, names);
    str = strtok( buf, seps );
    if( str == NULL ){
        printf( "\n You must determine the parameters names\n");
        return;
    }
    do {
       strcpy(str1, str);
       p = str1;
       for( i = 0; i < strlen(str1); i++ ) {
          *p = tolower(str1[i]);
          p++;
       }
       if ( strcmp( str1, "layout" ) == 0 ) {
           itmp = va_arg( ap, CBLAS_LAYOUT );
           if ( itmp == CblasRowMajor )
                printf("LAYOUT = CblasRowMajor  ");
           else if ( itmp == CblasColMajor )
                printf("LAYOUT = CblasColMajor  ");
       } else if ( strcmp( str1, "side" ) == 0 ) {
            itmp = va_arg( ap, CBLAS_SIDE);
            if ( itmp == CblasLeft )
               printf("SIDE = CblasLeft  ");
            else if ( itmp == CblasRight )
               printf("SIDE = CblasRight  ");
       } else if ( strcmp( str1, "uplo" ) == 0 ) {
            itmp = va_arg( ap, CBLAS_UPLO );
            if ( itmp == CblasUpper )
               printf("UPLO = CblasUpper  ");
            else if ( itmp == CblasLower )
               printf("UPLO = CblasLower  ");
       } else if ( strcmp( str1, "diag" ) == 0 ) {
            itmp = va_arg( ap, CBLAS_DIAG );
            if ( itmp == CblasUnit )
               printf("DIAG=CblasUnit  ");
            else if ( itmp == CblasNonUnit )
               printf("DIAG = CblasNonUnit  ");
       } else if ( ( strcmp( str1, "trans"  ) == 0 ) ||
                   ( strcmp( str1, "transa" ) == 0 ) ||
                   ( strcmp( str1, "transb" ) == 0 ) ) {
            itmp = va_arg( ap, CBLAS_TRANSPOSE );
            if ( itmp == CblasNoTrans )
               printf("%s = CblasNoTrans  ", str);
            else if ( itmp == CblasTrans )
               printf("%s = CblasTrans  ", str);
            else if ( itmp == CblasConjTrans )
               printf("%s = CblasConjTrans  ", str);
       } else if ( strcmp( str1, "offsetc" ) == 0 ) {
           itmp = va_arg( ap, CBLAS_OFFSET );
           if ( itmp == CblasFixOffset ) {
               printf("OFFSETC = CblasFixOffset ");
           } else if ( itmp == CblasRowOffset ) {
               printf("OFFSETC = CblasRowOffset ");
           } else if ( itmp == CblasColOffset ) {
               printf("OFFSETC = CblasColOffset ");
           }
       }
    } while ( (str = strtok( NULL, seps )) != NULL );
    va_end( ap );
    return;
} /* PrintParameters */

void PrintVectorI(MKL_INT n, MKL_INT *x, char *name)
{
      MKL_INT i;

      printf("\n       VECTOR %s\n         ", name);
      for (i = 0; i < n; i++)
          printf("    "INT_FORMAT, *(x+i));
      return;
} /* PrintVectorI */

void PrintVectorI16(MKL_INT n, MKL_INT16 *x, char *name)
{
      MKL_INT i;

      printf("\n       VECTOR %s\n         ", name);
      for (i = 0; i < n; i++)
          printf("    %d", *(x+i));
      return;
} /* PrintVectorI16 */

void PrintVectorI32(MKL_INT n, MKL_INT32 *x, char *name)
{
      MKL_INT i;

      printf("\n       VECTOR %s\n         ", name);
      for (i = 0; i < n; i++)
          printf("    %d", *(x+i));
      return;
} /* PrintVectorI32 */

void PrintVectorS(int flag, MKL_INT n, float *x,  MKL_INT incx, char *name)
{
      MKL_INT i;

      switch(flag) {
      case 0:
           printf("\n       VECTOR %s   INC%s=" INT_FORMAT "\n         ",
                    name, name, incx);
           break;
      case 1:
           printf("\n       VECTOR %s\n         ", name);
           break;
      default:
           break;
      } /* switch */
      for (i = 0; i < (1+(n-1)*ABS(incx)); i++)
          printf("%6.2f  ", *(x+i));
      return;
} /* PrintVectorS */

void PrintVectorD(int flag, MKL_INT n, double *x,  MKL_INT incx, char *name)
{
      MKL_INT i;

      switch(flag) {
      case 0:
           printf("\n       VECTOR %s   INC%s=" INT_FORMAT "\n         ",
                    name, name, incx);
           break;
      case 1:
           printf("\n       VECTOR %s\n         ", name);
           break;
      default:
           break;
      } /* switch */
      for (i = 0; i < (1+(n-1)*ABS(incx)); i++)
          printf("%8.3f  ", *(x+i));
      return;
} /* PrintVectorD */

void PrintVectorC(int flag, MKL_INT n, MKL_Complex8 *x,  MKL_INT incx, char *name)
{
      MKL_INT i;

      switch(flag) {
      case 0:
           printf("\n       VECTOR %s   INC%s=" INT_FORMAT "\n         ",
                    name, name, incx);
           break;
      case 1:
           printf("\n       VECTOR %s\n         ", name);
           break;
      default:
           break;
      } /* switch */
      for (i = 0; i < (1+(n-1)*ABS(incx)); i++)
          printf("(%6.2f,%6.2f)   ", (x+i)->real, (x+i)->imag);
      return;
} /* PrintVectorC */

void PrintVectorZ(int flag, MKL_INT n, MKL_Complex16 *x,  MKL_INT incx, char *name)
{
      MKL_INT i;

      switch(flag) {
      case 0:
           printf("\n       VECTOR %s   INC%s=" INT_FORMAT "\n         ",
                    name, name, incx);
               break;
      case 1:
           printf("\n       VECTOR %s\n         ", name);
               break;
      default:
               break;
      } /* switch */
      for (i = 0; i < (1+(n-1)*ABS(incx)); i++)
          printf("(%6.2f,%6.2f)   ", (x+i)->real, (x+i)->imag);
      return;
} /* PrintVectorZ */

void PrintArrayI8(CBLAS_LAYOUT*layout, int flag1, int flag2, MKL_INT *m, MKL_INT *n, MKL_INT8 *a,
                MKL_INT *lda, char *name)
{
        MKL_INT i, j;
        MKL_INT8 *addr;

        switch(flag1) {
        case 0:
            printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, *lda);
            break;
        case 1:
            printf("\n       ARRAY %s", name);
            break;
        default:
            break;
        } /* switch */
        if (*layout == CblasRowMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j < *n; j++)
                        printf("%d  ", (int) *(addr+j));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%d  ", (int) *(addr+j));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j <= i; j++)
                        printf("%d  ", (int) *(addr+j));
                } /* for */
            } /* if */
        } else if (*layout == CblasColMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j < *n; j++)
                        printf("%d  ", (int) *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%d  ", (int)*(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j <= i; j++)
                        printf("%d  ", (int) *(addr+j*(*lda)));
                } /* for */
            } /* if */
        } /* if */
        return;
} /* PrintArrayI8 */

void PrintArrayI16(CBLAS_LAYOUT*layout, int flag1, int flag2, MKL_INT *m, MKL_INT *n, MKL_INT16 *a,
                MKL_INT *lda, char *name)
{
        MKL_INT i, j;
        MKL_INT16 *addr;

        switch(flag1) {
        case 0:
            printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, *lda);
            break;
        case 1:
            printf("\n       ARRAY %s", name);
            break;
        default:
            break;
        } /* switch */
        if (*layout == CblasRowMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j < *n; j++)
                        printf("%hd  ", *(addr+j));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%hd  ", *(addr+j));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j <= i; j++)
                        printf("%hd  ", *(addr+j));
                } /* for */
            } /* if */
        } else if (*layout == CblasColMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j < *n; j++)
                        printf("%hd  ", *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%hd  ", *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j <= i; j++)
                        printf("%hd  ", *(addr+j*(*lda)));
                } /* for */
            } /* if */
        } /* if */
        return;
} /* PrintArrayI16 */

void PrintArrayI32(CBLAS_LAYOUT*layout, int flag1, int flag2, MKL_INT *m, MKL_INT *n, MKL_INT32 *a,
                MKL_INT *lda, char *name)
{
        MKL_INT i, j;
        MKL_INT32 *addr;

        switch(flag1) {
        case 0:
            printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, *lda);
            break;
        case 1:
            printf("\n       ARRAY %s", name);
            break;
        default:
            break;
        } /* switch */
        if (*layout == CblasRowMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j < *n; j++)
                        printf("%d  ", *(addr+j));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%d  ", *(addr+j));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j <= i; j++)
                        printf("%d  ", *(addr+j));
                } /* for */
            } /* if */
        } else if (*layout == CblasColMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j < *n; j++)
                        printf("%d  ", *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%d  ", *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j <= i; j++)
                        printf("%d  ", *(addr+j*(*lda)));
                } /* for */
            } /* if */
	} /* if */
	return;
} /* PrintArrayI32 */

void PrintArrayBF16(CBLAS_LAYOUT*layout, int flag1, int flag2, MKL_INT *m, MKL_INT *n, MKL_BF16 *a,
                MKL_INT *lda, char *name)
{
        MKL_INT i, j;
        MKL_BF16 *addr;

        switch(flag1) {
        case 0:
            printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, *lda);
            break;
        case 1:
            printf("\n       ARRAY %s", name);
            break;
        default:
            break;
        } /* switch */
        if (*layout == CblasRowMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j < *n; j++)
                        printf("%6.2f  ", b2f(*(addr+j)));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%6.2f  ", b2f(*(addr+j)));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j <= i; j++)
                        printf("%6.2f  ", b2f(*(addr+j)));
                } /* for */
            } /* if */
        } else if (*layout == CblasColMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j < *n; j++)
                        printf("%6.2f  ", b2f(*(addr+j*(*lda))));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%6.2f  ", b2f(*(addr+j*(*lda))));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j <= i; j++)
                        printf("%6.2f  ", b2f(*(addr+j*(*lda))));
                } /* for */
            } /* if */
        } /* if */
        return;
} /* PrintArrayBF16 */

void PrintArrayS(CBLAS_LAYOUT*layout, int flag1, int flag2, MKL_INT *m, MKL_INT *n, float *a,
                MKL_INT *lda, char *name)
{
        MKL_INT i, j;
        float   *addr;

        switch(flag1) {
        case 0:
            printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, *lda);
            break;
        case 1:
            printf("\n       ARRAY %s", name);
            break;
        default:
            break;
        } /* switch */
        if (*layout == CblasRowMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j < *n; j++)
                        printf("%6.2f  ", *(addr+j));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%6.2f  ", *(addr+j));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j <= i; j++)
                        printf("%6.2f  ", *(addr+j));
                } /* for */
            } /* if */
        } else if (*layout == CblasColMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j < *n; j++)
                        printf("%6.2f  ", *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j=0; j < i; j++)
                        printf("        ");
                    for (j = i; j < *n; j++)
                        printf("%6.2f  ", *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < (*m); i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j <= i; j++)
                        printf("%6.2f  ", *(addr+j*(*lda)));
                } /* for */
            } /* if */
        } /* if */
        return;
} /* PrintArrayS */

void PrintArrayD(CBLAS_LAYOUT*layout, int flag1, int flag2, MKL_INT *m, MKL_INT *n, double *a,
                MKL_INT *lda, char *name)
{
        MKL_INT i, j;
        double  *addr;

        switch(flag1) {
        case 0:
            printf("\n       ARRAY %s   LD%s=" INT_FORMAT,
                    name, name, *lda);
            break;
        case 1:
            printf("\n       ARRAY %s", name);
            break;
        default:
            break;
        } /* switch */
        if (*layout == CblasRowMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j < *n; j++)
                        printf("%8.3f  ", *(addr+j));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j=0; j < i; j++)
                        printf("          ");
                    for (j = i; j < *n; j++)
                        printf("%8.3f  ", *(addr+j));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j <= i; j++)
                        printf("%8.3f  ", *(addr+j));
                } /* for */
            } /* if */
        } else if (*layout == CblasColMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j < *n; j++)
                        printf("%8.3f  ", *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j=0; j < i; j++)
                        printf("          ");
                    for (j = i; j < *n; j++)
                        printf("%8.3f  ", *(addr+j*(*lda)));
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j <= i; j++)
                        printf("%8.3f  ", *(addr+j*(*lda)));
                } /* for */
            } /* if */
        } /* if */
        return;
} /* PrintArrayD */

void PrintArrayC(CBLAS_LAYOUT*layout, int flag1, int flag2, MKL_INT *m, MKL_INT *n, MKL_Complex8 *a,
                MKL_INT *lda, char *name)
{
        MKL_INT  i, j;
        MKL_Complex8 *addr;

        switch(flag1) {
        case 0:
            printf("\n       ARRAY %s   LD%s=" INT_FORMAT,
                    name, name, *lda);
            break;
        case 1:
            printf("\n       ARRAY %s", name);
            break;
        default:
            break;
        } /* switch */
        if (*layout == CblasRowMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j < *n; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j)->real, (addr+j)->imag);
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j=0; j < i; j++)
                        printf("                  ");
                    for (j = i; j < *n; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j)->real, (addr+j)->imag);
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j <= i; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j)->real, (addr+j)->imag);
                } /* for */
            } /* if */
        } else if(*layout == CblasColMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j < *n; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j*(*lda))->real,
                                (addr+j*(*lda))->imag);
                }  /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j=0; j < i; j++)
                        printf("                  ");
                    for (j = i; j < *n; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j*(*lda))->real,
                                (addr+j*(*lda))->imag);
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j <= i; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j*(*lda))->real,
                                (addr+j*(*lda))->imag);
                } /* for */
            } /* if */
        } /* if */
        return;
} /* PrintArrayC */

void PrintArrayZ(CBLAS_LAYOUT*layout, int flag1, int flag2,
                 MKL_INT *m, MKL_INT *n, MKL_Complex16 *a, MKL_INT *lda, char *name)
{
        MKL_INT   i, j;
        MKL_Complex16 *addr;

        if (flag1 == 0)
            printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, *lda);
        else
            printf("\n       ARRAY %s", name);

        if (*layout == CblasRowMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j < *n; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j)->real,
                               (addr+j)->imag);
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j=0; j < i; j++)
                        printf("                  ");
                    for (j = i; j < *n; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j)->real,
                               (addr+j)->imag);
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i*(*lda);
                    for (j = 0; j <= i; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j)->real,
                                (addr+j)->imag);
                } /* for */
            } /* if */
        } else if (*layout == CblasColMajor) {
            if (flag2 == GENERAL_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j < *n; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j*(*lda))->real,
                                (addr+j*(*lda))->imag);
                } /* for */
            } else if (flag2 == UPPER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j=0; j < i; j++)
                        printf("                  ");
                    for (j = i; j < *n; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j*(*lda))->real,
                                (addr+j*(*lda))->imag);
                } /* for */
            } else if (flag2 == LOWER_MATRIX) {
                for (i = 0; i < *m; i++) {
                    printf("\n         ");
                    addr = a + i;
                    for (j = 0; j <= i; j++)
                        printf("(%6.2f,%6.2f)   ", (addr+j*(*lda))->real,
                                (addr+j*(*lda))->imag);
                } /* for */
            } /* if */
        } /* if */
        return;
} /* PrintArrayZ */

void PrintBandArrayS(CBLAS_LAYOUT*layout, int flag,
                     MKL_INT kl, MKL_INT ku, MKL_INT m, MKL_INT n,
                     float *a, MKL_INT lda, char *name)
{
      MKL_INT    i, j, l;
      float *addr, *addr1;
      MKL_INT    kl1, ku1, i_start, j_start, j_end, kl_rows, ku_rows;

      if (flag == 0)
          printf("\n       ARRAY %s   LD%s=" INT_FORMAT "  KL=" INT_FORMAT "  KU=" INT_FORMAT,
                   name, name, lda, kl, ku);
      else
          printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, lda);

      if (*layout == CblasRowMajor) {
         for( i = 0; i < MIN( m, n ); i++ ) {
            printf("\n    ");
            addr = a + i*lda;
            kl1  = ( i - kl <= 0 ) ? i : kl;
            ku1  = ( i + ku >= n ) ? MAX(0,n-i-1) : ku;
            j_start = kl - kl1;
            j_end = j_start + kl1 + ku1;
            addr1 = addr + j_start;
            for( l=0; l < j_start; l++ )
                printf("        ");
            for( l=0; l < j_end-j_start+1; l++ )
                printf("%6.2f  ", *(addr1+l));
         }
         for( i = MIN( m, n ); i < MIN( m, MIN( m, n ) + kl); i++ ) {
            printf("\n    ");
            addr = a + i*lda;
            kl1  = n - i + kl;
            j_start = ( kl > n ) ? kl - n : n - kl;
            j_end = j_start + kl1 - 1;
            addr1 = addr + j_start;
            for( l=0; l < j_start; l++ )
                printf("        ");
            for( l=0; l < j_end-j_start+1; l++ )
                printf("%6.2f  ", *(addr1+l));
         }
      } else if (*layout == CblasColMajor) {
         i_start = (ku > n ) ? ku - n + 1 : 0;
         ku_rows = (ku > n ) ? n - 1 : ku;
         j_start = ku_rows;
         for( i = 0; i < ku_rows; i++ ) {
              printf("\n      ");
              j = j_start*lda; addr1 = a + i + i_start;
              for( l=0; l < j_start; l++ )
                 printf("        ");
              for( l=0; l < n-ku_rows; l++ )
                 printf("%6.2f  ", *(addr1+j+l*lda) );
              j_start--;
         }

         j_end = MIN(m,n);
         addr1 = a + ku;
         printf("\n      ");
         for( l=0; l < j_end; l++ )
            printf("%6.2f  ", *(addr1+l*lda) );

         kl_rows = ( kl <= m-1 ) ?  kl : m - 1;
         for ( i = 1; i < kl_rows+1; i++ ) {
              printf("\n      ");
              kl1 = ( i+j_end <= m ) ? j_end : m - i;
              addr1 = a + ku + i;
              for( l=0; l < kl1; l++ )
                 printf("%6.2f  ", *(addr1+l*lda) );
         }
      } /* if */
      return;
} /* PrintBandArrayS */

void PrintBandArrayD(CBLAS_LAYOUT*layout, int flag,
                     MKL_INT kl, MKL_INT ku, MKL_INT m, MKL_INT n,
                     double *a, MKL_INT lda, char *name)
{
      MKL_INT i, j, l;
      double *addr, *addr1;
      MKL_INT kl1, ku1, i_start, j_start, j_end, kl_rows, ku_rows;

      if (flag == 0)
          printf("\n       ARRAY %s   LD%s=" INT_FORMAT "  KL=" INT_FORMAT "  KU=" INT_FORMAT,
                   name, name, lda, kl, ku);
      else
          printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, lda);

      if (*layout == CblasRowMajor) {
         for( i = 0; i < MIN( m, n ); i++ ) {
            printf("\n    ");
            addr = a + i*lda;
            kl1  = ( i - kl <= 0 ) ? i : kl;
            ku1  = ( i + ku >= n ) ? MAX(0,n-i-1) : ku;
            j_start = kl - kl1;
            j_end = j_start + kl1 + ku1;
            addr1 = addr + j_start;
            for( l=0; l < j_start; l++ )
                printf("        ");
            for( l=0; l < j_end-j_start+1; l++ )
                printf("%6.2f  ", *(addr1+l));
         }
         for( i = MIN( m, n ); i < MIN( m, MIN( m, n ) + kl); i++ ) {
            printf("\n    ");
            addr = a + i*lda;
            kl1  = n - i + kl;
            j_start = ( kl > n ) ? kl - n : n - kl;
            j_end = j_start + kl1 - 1;
            addr1 = addr + j_start;
            for( l=0; l < j_start; l++ )
                printf("        ");
            for( l=0; l < j_end-j_start+1; l++ )
                printf("%6.2f  ", *(addr1+l));
         }
      } else if (*layout == CblasColMajor) {
         i_start = (ku > n ) ? ku - n + 1 : 0;
         ku_rows = (ku > n ) ? n - 1 : ku;
         j_start = ku_rows;
         for( i = 0; i < ku_rows; i++ ) {
              printf("\n      ");
              j = j_start*lda; addr1 = a + i + i_start;
              for( l=0; l < j_start; l++ )
                 printf("        ");
              for( l=0; l < n-ku_rows; l++ )
                 printf("%6.2f  ", *(addr1+j+l*lda) );
              j_start--;
         }

         j_end = MIN(m,n);
         addr1 = a + ku;
         printf("\n      ");
         for( l=0; l < j_end; l++ )
            printf("%6.2f  ", *(addr1+l*lda) );

         kl_rows = ( kl <= m-1 ) ?  kl : m - 1;
         for ( i = 1; i < kl_rows+1; i++ ) {
              printf("\n      ");
              kl1 = ( i+j_end <= m ) ? j_end : m - i;
              addr1 = a + ku + i;
              for( l=0; l < kl1; l++ )
                 printf("%6.2f  ", *(addr1+l*lda) );
         }
      } /* if */
      return;
} /* PrintBandArrayD */

void PrintBandArrayC(CBLAS_LAYOUT*layout, int flag,
                     MKL_INT kl, MKL_INT ku, MKL_INT m, MKL_INT n,
                     MKL_Complex8 *a, MKL_INT lda, char *name)
{
      MKL_INT       i, j, l;
      MKL_Complex8 *addr, *addr1;
      MKL_INT       kl1, ku1, j_start, j_end, i_start, kl_rows, ku_rows;

      if (flag == 0)
          printf("\n       ARRAY %s   LD%s=" INT_FORMAT "  KL=" INT_FORMAT "  KU=" INT_FORMAT,
                   name, name, lda, kl, ku);
      else
          printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, lda);

      if (*layout == CblasRowMajor) {
         for( i = 0; i < MIN( m, n ); i++ ) {
            printf("\n    ");
            addr = a + i*lda;
            kl1  = ( i - kl <= 0 ) ? i : kl;
            ku1  = ( i + ku >= n ) ? MAX(0,n-i-1) : ku;
            j_start = kl - kl1;
            j_end = j_start + kl1 + ku1;
            addr1 = addr + j_start;
            for( l=0; l < j_start; l++ )
                printf("                 ");
            for( l=0; l < j_end-j_start+1; l++ )
                printf("(%6.2f,%6.2f)  ", (addr1+l)->real, (addr1+l)->imag);
         }
         for( i = MIN( m, n ); i < MIN( m, MIN( m, n ) + kl); i++ ) {
            printf("\n    ");
            addr = a + i*lda;
            kl1  = n - i + kl;
            j_start = ( kl > n ) ? kl - n : n - kl;
            j_end = j_start + kl1 - 1;
            addr1 = addr + j_start;
            for( l=0; l < j_start; l++ )
                printf("                 ");
            for( l=0; l < j_end-j_start+1; l++ )
                printf("(%6.2f,%6.2f)  ", (addr1+l)->real, (addr1+l)->imag);
         }
      } else if (*layout == CblasColMajor) {
         i_start = (ku > n ) ? ku - n + 1 : 0;
         ku_rows = (ku > n ) ? n - 1 : ku;
         j_start = ku_rows;
         for( i = 0; i < ku_rows; i++ ) {
              printf("\n               ");
              j = j_start*lda; addr1 = a + i + i_start;
              for( l=0; l < j_start; l++ )
                 printf("                 ");
              for( l=0; l < n-ku_rows; l++ )
                 printf("(%6.2f,%6.2f)  ",
                         (addr1+j+l*lda)->real, (addr1+j+l*lda)->imag );
              j_start--;
         }

         j_end = MIN(m,n);
         addr1 = a + ku;
         printf("\n               ");
         for( l=0; l < j_end; l++ )
            printf("(%6.2f,%6.2f)  ",
                    (addr1+l*lda)->real, (addr1+l*lda)->imag );

         kl_rows = ( kl <= m-1 ) ?  kl : m - 1;
         for ( i = 1; i < kl_rows+1; i++ ) {
              printf("\n               ");
              kl1 = ( i+j_end <= m ) ? j_end : m - i;
              addr1 = a + ku + i;
              for( l=0; l < kl1; l++ )
                 printf("(%6.2f,%6.2f)  ",
                         (addr1+l*lda)->real, (addr1+l*lda)->imag );
         }
      } /* if */
      return;
} /* PrintBandArrayC */

void PrintBandArrayZ(CBLAS_LAYOUT*layout, int flag,
                     MKL_INT kl, MKL_INT ku, MKL_INT m, MKL_INT n,
                     MKL_Complex16 *a, MKL_INT lda, char *name)
{
      MKL_INT        i, j, l;
      MKL_Complex16 *addr, *addr1;
      MKL_INT        kl1, ku1, j_start, j_end, i_start, kl_rows, ku_rows;

      if (flag == 0)
          printf("\n       ARRAY %s   LD%s=" INT_FORMAT "  KL=" INT_FORMAT "  KU=" INT_FORMAT,
                   name, name, lda, kl, ku);
      else
          printf("\n       ARRAY %s   LD%s=" INT_FORMAT, name, name, lda);

      if (*layout == CblasRowMajor) {
         for( i = 0; i < MIN( m, n ); i++ ) {
            printf("\n    ");
            addr = a + i*lda;
            kl1  = ( i - kl <= 0 ) ? i : kl;
            ku1  = ( i + ku >= n ) ? MAX(0,n-i-1) : ku;
            j_start = kl - kl1;
            j_end = j_start + kl1 + ku1;
            addr1 = addr + j_start;
            for( l=0; l < j_start; l++ )
                printf("                 ");
            for( l=0; l < j_end-j_start+1; l++ )
                printf("(%6.2f,%6.2f)  ", (addr1+l)->real, (addr1+l)->imag);
         }
         for( i = MIN( m, n ); i < MIN( m, MIN( m, n ) + kl); i++ ) {
            printf("\n    ");
            addr = a + i*lda;
            kl1  = n - i + kl;
            j_start = ( kl > n ) ? kl - n : n - kl;
            j_end = j_start + kl1 - 1;
            addr1 = addr + j_start;
            for( l=0; l < j_start; l++ )
                printf("                 ");
            for( l=0; l < j_end-j_start+1; l++ )
                printf("(%6.2f,%6.2f)  ", (addr1+l)->real, (addr1+l)->imag);
         }
      } else if (*layout == CblasColMajor) {
         i_start = (ku > n ) ? ku - n + 1 : 0;
         ku_rows = (ku > n ) ? n - 1 : ku;
         j_start = ku_rows;
         for( i = 0; i < ku_rows; i++ ) {
              printf("\n      ");
              j = j_start*lda; addr1 = a + i + i_start;
              for( l=0; l < j_start; l++ )
                 printf("                 ");
              for( l=0; l < n-ku_rows; l++ )
                 printf("(%6.2f,%6.2f)  ",
                         (addr1+j+l*lda)->real, (addr1+j+l*lda)->imag );
              j_start--;
         }

         j_end = MIN(m,n);
         addr1 = a + ku;
          printf("\n      ");
         for( l=0; l < j_end; l++ )
            printf("(%6.2f,%6.2f)  ",
                    (addr1+l*lda)->real, (addr1+l*lda)->imag );

         kl_rows = ( kl <= m-1 ) ?  kl : m - 1;
         for ( i = 1; i < kl_rows+1; i++ ) {
              printf("\n      ");
              kl1 = ( i+j_end <= m ) ? j_end : m - i;
              addr1 = a + ku + i;
              for( l=0; l < kl1; l++ )
                 printf("(%6.2f,%6.2f)  ",
                         (addr1+l*lda)->real, (addr1+l*lda)->imag );
         }
      } /* if */
      return;
} /* PrintBandArrayZ */

