/*******************************************************************************
* Copyright 2006-2020 Intel Corporation.
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
! Content:
!          Intel(R) Math Kernel Library (Intel(R) MKL) DFTI implementation
!          through FFTW interface (via wrappers) example program (C-interface)
!
! Complex-to-complex 1D transform for double/single precision data not inplace.
!
! Configuration parameters for Intel MKL DFTI:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX     (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE/DFTI_SINGLE (obligatory)
!           DFTI_DIMENSION      = 1                (obligatory)
!           DFTI_LENGTHS        = n                (obligatory)
!           DFTI_PLACEMENT      = DFTI_NOT_INPLACE (default=DFTI_INPLACE)
!           DFTI_FORWARD_SCALE  = 1.0              (default)
!           DFTI_BACKWARD_SCALE = 1.0/n            (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.h interface file
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fftw.h"
#include "mkl_dfti_examples.h"

int main(void)  /* COMPLEX_1D_EX4 */
{
    /*
    **  DFT input parameters
    */
    int n = 4;
    int multiple = 3;
    int istride = 2;
    int ostride = 2;
    int idist = 8;
    int odist = 10;

    int lda, lda_out;

    fftw_plan my_plan;

    fftw_complex* x_in;
    fftw_complex* x_out;
    fftw_complex* x_exp;

    TYPE_PRECISION Scale;
    TYPE_PRECISION maxerr;
    TYPE_PRECISION eps = EPS;

    if (LEGEND_PRINT) {
        printf("\n\n COMPLEX_1D_EX4\n");
        printf(" Forward-Backward 1D complex transform for double/single precision data\n\n");
        printf(" Configuration parameters:\n\n");
        printf(" DFTI_FORWARD_DOMAIN = DFTI_COMPLEX\n");
        printf(" DFTI_PRECISION      = DFTI_DOUBLE/DFTI_SINGLE\n");
        printf(" DFTI_DIMENSION      = 1\n");
        printf(" DFTI_LENGTHS        = %d\n", n);
        printf(" DFTI_PLACEMENT      = DFTI_NOT_INPLACE\n");
        printf(" DFTI_FORWARD_SCALE  = 1.0\n");
        printf(" DFTI_BACKWARD_SCALE = 1.0/n\n\n");
    }

    /*
    **  Allocate array for input/output and expected data
    */
    lda = n*istride*multiple;
    if (idist > n*istride) lda = idist*multiple;
    lda_out = n*ostride*multiple;
    if (odist > n*ostride) lda_out = odist*multiple;

    x_in = (fftw_complex*)fftw_malloc(2*lda*sizeof(TYPE_PRECISION));
    x_out = (fftw_complex*)fftw_malloc(2*lda_out*sizeof(TYPE_PRECISION));
    x_exp = (fftw_complex*)fftw_malloc(2*lda*sizeof(TYPE_PRECISION));

    /*
    **  Initialize x_in and copy to expected x_exp
    */
#ifdef MKL_DOUBLE
    init_input_and_expected_vectors_z(x_in, x_exp, lda);
#else
    init_input_and_expected_vectors_c(x_in, x_exp, lda);
#endif

    /*
    **  Create FFTW plan for 1D double/single precision forward transform
    */
    my_plan = fftw_create_plan(n, FFTW_FORWARD, FFTW_ESTIMATE);
    DIE_UNLESS(my_plan);

    /*
    **  Compute DFT
    */
    fftw(my_plan, multiple, x_in, istride, idist, x_out, ostride, odist);

#ifdef DEBUG
    {
        int i;
        for (i=0; i<lda; i++) {
            printf(" forward: x_out[%d].re = %f, x_out[%d].im = %f\n", i, x_out[i].re, i, x_out[i].im);
        }
    }
#endif

    /*
    **  Destroy FFTW plan
    */
    fftw_destroy_plan(my_plan);

    /*
    **  Set Scale number for Backward transform
    */
    Scale = 1.0/(TYPE_PRECISION)n;

    /*
    **  Create FFTW plan for 1D double/single precision backward transform
    */
    my_plan = fftw_create_plan(n, FFTW_BACKWARD, FFTW_ESTIMATE);
    DIE_UNLESS(my_plan);

    /*
    **  Compute DFT
    */
    fftw(my_plan, multiple, x_out, ostride, odist, x_in, istride, idist);

#ifdef DEBUG
    {
        int i;
        for (i=0; i<lda; i++) {
            printf(" backward: x_in[%d].re = %f, x_in[%d].im = %f\n", i, x_in[i].re, i, x_in[i].im);
        }
    }
#endif

    /*
    **  Destroy FFTW plan
    */
    fftw_destroy_plan(my_plan);

    /*
    **  Result scaling
    */
    scaling_multiple(x_in, Scale, n, multiple, istride, idist);

    /*
    **  Check result
    */
#ifdef MKL_DOUBLE
    maxerr = check_result_z(x_in, x_exp, lda);
#else
    maxerr = check_result_c(x_in, x_exp, lda);
#endif

    /*
    **  Free arrays for input/output and expected data
    */
    fftw_free(x_in);
    fftw_free(x_out);
    fftw_free(x_exp);

    if (ACCURACY_PRINT)
        printf("\n Accuracy = %g\n\n", maxerr);

    if (maxerr < eps) {
        printf(" TEST PASSED\n");
    } else {
        printf(" TEST FAILED\n");
        return 1;
    }

    printf(" END OF TEST\n");

    return 0;
}
