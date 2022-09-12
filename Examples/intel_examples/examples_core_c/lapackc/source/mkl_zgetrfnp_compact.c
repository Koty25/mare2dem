/*******************************************************************************
* Copyright 2017-2020 Intel Corporation.
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
   MKL_ZGETRFNP_COMPACT Example.
   ==============

   Program computes the LU factorization, without pivoting, of a set of matrices
   A1, ... , Anmat using compact routines. 
   For this example, nmat = 512; however, only the results and input for the
   first matrix A1 are shown.

   Description.
   ============

   The routine computes the LU factorization, without pivoting, of a set of
   general, m x n matrices, that have been stored in compact format. 
*/
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#define LAYOUT    MKL_COL_MAJOR
#define M                     5
#define N                     4
#define NMAT                512

/* Auxiliary routines prototypes */
extern void print_matrix(MKL_LAYOUT layout, MKL_INT m, MKL_INT n, MKL_Complex16* a,
                         MKL_INT lda);

int main() {
    MKL_INT i, j;

    MKL_LAYOUT layout = LAYOUT;
    MKL_INT col_major = (layout == MKL_COL_MAJOR);
    MKL_INT m = M;
    MKL_INT n = N;
    MKL_INT lda = col_major ? m : n;
    MKL_INT info;
    MKL_COMPACT_PACK format;
    MKL_INT nmat = NMAT;

    /* For setting up compact arrays */
    MKL_INT a_buffer_size;
    MKL_INT ldap = lda;
    MKL_INT sdap = (col_major ? n : m);
    MKL_INT a_size = lda * (col_major ? n : m);

    /* Set up standard arrays in P2P (pointer-to-pointer) format */
    MKL_INT na = lda * (col_major ? n : m) * nmat;
    MKL_Complex16 *a = (MKL_Complex16 *)mkl_malloc(na * nmat * sizeof(MKL_Complex16), 128);
    MKL_Complex16 *a_array[NMAT];

    double *a_compact;

    /* For random generation of matrices */
    MKL_INT idist = 1;
    MKL_INT iseed[] = { 0, 1, 2, 3 };
    double diag_offset = (double)m;

    /* Random generation of matrices */
    zlarnv(&idist, iseed, &na, a);

    for (i = 0; i < nmat; i++) {
        /* Make matrix diagonal dominant to avoid accuracy issues
                 in the non-pivoted LU factorization */
        for (j = 0; j < n; j++) {
            a[i * a_size + j + j * lda].real += diag_offset;
            a[i * a_size + j + j * lda].imag = 0.0f;
        }
        a_array[i] = &a[i * a_size];
    }

    /* Print input data */
    printf("   LAYOUT=%s, M=%d, N=%d, NMAT=%d\n",
            (col_major ? "COL-MAJOR" : "ROW-MAJOR"), m, n, nmat);
    printf("\n   Matrix A1:\n");
    print_matrix(layout, m, n, a_array[0], lda);

    /* Set up Compact arrays */
    format = mkl_get_format_compact();

    a_buffer_size = mkl_zget_size_compact(ldap, sdap, format, nmat);

    a_compact = (double *)mkl_malloc(a_buffer_size, 128);

    /* Pack from P2P to Compact format */
    mkl_zgepack_compact(layout, m, n, a_array, lda, a_compact, ldap, format, nmat);

    /* Perform Compact LU Factorization */
    mkl_zgetrfnp_compact(layout, m, n, a_compact, ldap, &info, format, nmat);

    /* Unpack from Compact to P2P format */
    mkl_zgeunpack_compact(layout, m, n, a_array, lda, a_compact, ldap, format, nmat);

    /* Print output data */
    printf("\n\n   Matrix A1 after factorization:\n");
    print_matrix(layout, m, n, a_array[0], lda);

    mkl_free(a_compact);
    mkl_free(a);

    return 0;
}

/* Auxiliary routine: printing a matrix */
void print_matrix(MKL_LAYOUT layout, MKL_INT m, MKL_INT n, MKL_Complex16* a, MKL_INT lda) {
    MKL_INT i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (layout == MKL_COL_MAJOR) {
                printf(" (%6.2f,%6.2f)", a[i + j * lda].real, a[i + j * lda].imag);
            } else {
                printf(" (%6.2f,%6.2f)", a[j + i * lda].real, a[j + i * lda].imag);
            }
        }
        printf("\n");
    }
}
