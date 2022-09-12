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
   MKL_DGETRINP_COMPACT Example.
   ==============

   Program computes the inverse, without pivoting, of a set of matrices A1,
   ... , Anmat, using compact routines. The LU-factorization of each Ai,
   1 <= i <= nmat, is found, and then the inverse of each Ai is calculated from
   the LU factorization. 
   For this example, nmat = 512; however, only the results and input for the 
   first matrix A1 are shown.

   Description.
   ============

   The routine computes the inverse, of a set of LU-factorized, without
   pivoting, general matrices, which are stored in the compact format.
*/
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#define LAYOUT    MKL_COL_MAJOR
#define N                     5
#define NMAT                512

/* Auxiliary routines prototypes */
extern void print_matrix(MKL_LAYOUT layout, MKL_INT m, MKL_INT n, double* a,
                         MKL_INT lda);

int main() {
    MKL_INT i, j;

    MKL_LAYOUT layout = LAYOUT;
    MKL_INT col_major = (layout == MKL_COL_MAJOR);
    MKL_INT n = N;
    MKL_INT lda = n;
    MKL_INT info;
    MKL_COMPACT_PACK format;
    MKL_INT nmat = NMAT;

    /* For setting up compact arrays */
    MKL_INT a_buffer_size;
    MKL_INT ldap = lda;
    MKL_INT sdap = n;
    MKL_INT a_size = lda * n;

    /* Set up standard arrays in P2P (pointer-to-pointer) format */
    MKL_INT na = lda * n * nmat;
    double *a = (double *)mkl_malloc(na * nmat * sizeof(double), 128);
    double *a_array[NMAT];

    double *a_compact;

    /* For random generation of matrices */
    MKL_INT idist = 1;
    MKL_INT iseed[] = { 0, 1, 2, 3 };
    double diag_offset = (double)n;

    /* For workspace calculation */
    MKL_INT lwork = -1;
    double work_query[1];
    double *work;

    /* Random generation of matrices */
    dlarnv(&idist, iseed, &na, a);

    for (i = 0; i < nmat; i++) {
        /* Make matrix diagonal dominant to avoid accuracy issues
                 in the non-pivoted LU factorization */
        for (j = 0; j < n; j++) {
            a[i * a_size + j + j * lda] += diag_offset;
        }
        a_array[i] = &a[i * a_size];
    }

    /* Print input data */
    printf("   LAYOUT=%s, N=%d, NMAT=%d\n",
            (col_major ? "COL-MAJOR" : "ROW-MAJOR"), n, nmat);
    printf("\n   Matrix A1:\n");
    print_matrix(layout, n, n, a_array[0], lda);

    /* Set up Compact arrays */
    format = mkl_get_format_compact();

    a_buffer_size = mkl_dget_size_compact(ldap, sdap, format, nmat);

    a_compact = (double *)mkl_malloc(a_buffer_size, 128);

    /* Pack from P2P to Compact format */
    mkl_dgepack_compact(layout, n, n, a_array, lda, a_compact, ldap, format, nmat);

    /* Workspace Query */
    mkl_dgetrinp_compact(layout, n, a_compact, ldap, work_query, lwork, &info, format, nmat);
    lwork = (MKL_INT)work_query[0];
    work = (double *)mkl_malloc(sizeof(double) * lwork, 128);

    /* Perform Compact LU Factorization */
    mkl_dgetrfnp_compact(layout, n, n, a_compact, ldap, &info, format, nmat);

    /* Perform Compact Inverse Calculation */
    mkl_dgetrinp_compact(layout, n, a_compact, ldap, work, lwork, &info, format, nmat);

    /* Unpack from Compact to P2P format */
    mkl_dgeunpack_compact(layout, n, n, a_array, lda, a_compact, ldap, format, nmat);

    /* Print output data */
    printf("\n\n   Inverse(A1):\n");
    print_matrix(layout, n, n, a_array[0], lda);

    mkl_free(a);
    mkl_free(a_compact);
    mkl_free(work);

    return 0;
}

/* Auxiliary routine: printing a matrix */
void print_matrix(MKL_LAYOUT layout, MKL_INT m, MKL_INT n, double* a, MKL_INT lda) {
    MKL_INT i, j;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            if (layout == MKL_COL_MAJOR) {
                printf(" %6.2f", a[i + j * lda]);
            } else {
                printf(" %6.2f", a[j + i * lda]);
            }
        }
        printf("\n");
    }
}
