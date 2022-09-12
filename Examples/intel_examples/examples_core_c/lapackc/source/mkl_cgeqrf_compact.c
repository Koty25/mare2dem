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
   MKL_CGEQRF_COMPACT Example.
   ==============

   Program computes the QR factorization of a set of matrices A1, ... , Anmat,
   using compact routines.
   For this example, nmat = 512; however, only the results and input for the
   first matrix A1 are shown.

   Description.
   ============

   The routine computes the QR factorization, of a set of general,
   m x n matrices, that have been stored in compact format. 
*/
#include <stdio.h>
#include <stdlib.h>
#include "mkl.h"

#define LAYOUT    MKL_COL_MAJOR
#define M                     5
#define N                     4
#define NMAT                512

/* Auxiliary routines prototypes */
extern void print_matrix(MKL_LAYOUT layout, MKL_INT m, MKL_INT n, MKL_Complex8* a,
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
    MKL_INT minmn = (m < n) ? m : n;

    /* For setting up compact arrays */
    MKL_INT a_buffer_size;
    MKL_INT tau_buffer_size;
    MKL_INT ldap = lda;
    MKL_INT sdap = col_major ? n : m;

    MKL_INT na = lda * (col_major ? n : m) * nmat;
    MKL_Complex8 *a = (MKL_Complex8 *)mkl_malloc(na * nmat * sizeof(MKL_Complex8), 128);
    MKL_Complex8 *tau = (MKL_Complex8 *)mkl_malloc(minmn * nmat * sizeof(MKL_Complex8), 128);
    MKL_INT a_size = lda * (col_major ? n : m);

    /* Set up standard arrays in P2P (pointer-to-pointer) format */
    MKL_Complex8 *a_array[NMAT];
    MKL_Complex8 *tau_array[NMAT];

    float *a_compact;
    float *tau_compact;

    float work_query[1];
    float *work;
    MKL_INT lwork = -1;

    /* For random generation of matrices */
    MKL_INT idist = 1;
    MKL_INT iseed[] = { 0, 1, 2, 3 };

    /* Random generation of matrices */
    clarnv(&idist, iseed, &na, a);

    for (i = 0; i < nmat; i++) {
        a_array[i] = &a[i * a_size];
        tau_array[i] = &tau[i * minmn];
    }

    /* Print input data */
    printf("   LAYOUT=%s, M=%d, N=%d, NMAT=%d\n",
            (col_major ? "COL-MAJOR" : "ROW-MAJOR"), m, n, nmat);
    printf("\n   Matrix A1:\n");
    print_matrix(layout, m, n, a_array[0], lda);

    /* Set up Compact arrays */
    format = mkl_get_format_compact();

    a_buffer_size = mkl_cget_size_compact(ldap, sdap, format, nmat);
    tau_buffer_size = mkl_cget_size_compact(minmn, 1, format, nmat);

    a_compact = (float *)mkl_malloc(a_buffer_size, 128);
    tau_compact = (float *)mkl_malloc(tau_buffer_size, 128);

    /* Pack A matrix from P2P to Compact format */
    mkl_cgepack_compact(layout, m, n, a_array, lda, a_compact, ldap, format, nmat);

    /* Workspace query */
    mkl_cgeqrf_compact(layout, m, n, a_compact, ldap, tau_compact, work_query, lwork, &info, format, nmat);

    lwork = (MKL_INT)work_query[0];
    work = (float *)mkl_malloc(lwork * sizeof(float), 128);

    /* Perform Compact QR Factorization */
    mkl_cgeqrf_compact(layout, m, n, a_compact, ldap, tau_compact, work, lwork, &info, format, nmat);

    /* Unpack A from Compact to P2P format */
    mkl_cgeunpack_compact(layout, m, n, a_array, lda, a_compact, ldap, format, nmat);

    /* Unpack tau from Compact to P2P format */
    mkl_cgeunpack_compact(MKL_COL_MAJOR, minmn, 1, tau_array, minmn, tau_compact, minmn, format, nmat);

    /* Print output data */
    printf("\n\n   Matrix A1 after factorization:\n");
    print_matrix(layout, m, n, a_array[0], lda);
    printf("\n   Array tau1 after factorization:\n");
    print_matrix(MKL_COL_MAJOR, 1, minmn, tau_array[0], 1);

    mkl_free(a);
    mkl_free(tau);
    mkl_free(a_compact);
    mkl_free(tau_compact);
    mkl_free(work);

    return 0;
}

/* Auxiliary routine: printing a matrix */
void print_matrix(MKL_LAYOUT layout, MKL_INT m, MKL_INT n, MKL_Complex8* a, MKL_INT lda) {
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
