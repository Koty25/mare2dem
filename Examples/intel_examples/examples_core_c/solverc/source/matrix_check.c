/*******************************************************************************
* Copyright 2004-2020 Intel Corporation.
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
*   Content : Intel(R) Math Kernel Library (Intel(R) MKL) PARDISO C example
*
********************************************************************************
*/
/* -------------------------------------------------------------------- */
/* Example program to show the use of the "matrix_check" routine */
/* on symmetric linear systems */
/* -------------------------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mkl_types.h"
#include "mkl_sparse_handle.h"

int main (void)
{
    /* Matrix data. */
    MKL_INT n = 8;
    MKL_INT ia[9] = { 1, 5, 8, 10, 12, 15, 17, 18, 19};
    MKL_INT ja[18] =
      { 1,    3,       6, 7,
           2, 3,    5,
              3,             8,
                 4,       7,
                 4,    6, 7,
                       6,    8,
                          7,
                             8
      };
    double a[18] =
      { 7.0,      1.0,           2.0, 7.0,
            -4.0, 8.0,      2.0,
                  1.0,                     5.0,
                       7.0,           9.0,
                       5.0,      1.0, 5.0,
                                -1.0,      5.0,
                                     11.0,
                                           5.0
      };

    sparse_checker_error_values check_err_val;
    sparse_struct pt;
    int error = 0;


    sparse_matrix_checker_init(&pt);
    pt.n = n;
    pt.csr_ia = ia;
    pt.csr_ja = ja;
    pt.indexing         = MKL_ONE_BASED;
    pt.matrix_structure = MKL_UPPER_TRIANGULAR;
    pt.print_style      = MKL_C_STYLE;
    pt.message_level    = MKL_PRINT;

    check_err_val = sparse_matrix_checker(&pt);

    printf("Matrix check details: (%d, %d, %d)\n", pt.check_result[0], pt.check_result[1], pt.check_result[2]);

    if ( check_err_val == MKL_SPARSE_CHECKER_NONTRIANGULAR) {
        printf("Matrix check result: MKL_SPARSE_CHECKER_NONTRIANGULAR\n");
        error = 0;
    }
    else {
        if ( check_err_val == MKL_SPARSE_CHECKER_SUCCESS) { printf("Matrix check result: MKL_SPARSE_CHECKER_SUCCESS\n"); }
        if ( check_err_val == MKL_SPARSE_CHECKER_NON_MONOTONIC) { printf("Matrix check result: MKL_SPARSE_CHECKER_NON_MONOTONIC\n"); }
        if ( check_err_val == MKL_SPARSE_CHECKER_OUT_OF_RANGE) { printf("Matrix check result: MKL_SPARSE_CHECKER_OUT_OF_RANGE\n"); }
        if ( check_err_val == MKL_SPARSE_CHECKER_NONORDERED) { printf("Matrix check result: MKL_SPARSE_CHECKER_NONORDERED\n"); }
        error = 1;
    }

    return error;
}
