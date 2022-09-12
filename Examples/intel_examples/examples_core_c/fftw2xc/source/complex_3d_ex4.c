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
! Forward-Backward 3D complex transform for double/single precision data not inplace.
!
! Configuration parameters:
!           DFTI_FORWARD_DOMAIN = DFTI_COMPLEX     (obligatory)
!           DFTI_PRECISION      = DFTI_DOUBLE/DFTI_SINGLE (obligatory)
!           DFTI_DIMENSION      = 3                (obligatory)
!           DFTI_LENGTHS        = {m,n,k}          (obligatory)
!           DFTI_PLACEMENT      = DFTI_NOT_INPLACE (default=DFTI_INPLACE)
!           DFTI_INPUT_STRIDES  = is               (default={0,n*k,k,1})
!           DFTI_OUTPUT_STRIDES = os               (default={0,n*k,k,1})
!           DFTI_FORWARD_SCALE  = 1.0              (default)
!           DFTI_BACKWARD_SCALE = 1.0/(m*n*k)      (default=1.0)
!
! Other default configuration parameters are in the mkl_dfti.h interface file
!****************************************************************************/

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "fftw.h"
#include "mkl_dfti_examples.h"

int main(void)  /* COMPLEX_3D_EX4 */
{
    /*
    **  DFT input parameters
    */
    int m = 4;
    int n = 5;
    int k = 6;
    int multiple = 4;
    int istride = 2;
    int idist = 800;
    int ostride = 2;
    int odist = 750;
    int lda, lda_out;

    fftwnd_plan my_plan;
    fftw_complex* x_in;
    fftw_complex* x_out;
    fftw_complex* x_exp;

    TYPE_PRECISION Scale;
    TYPE_PRECISION maxerr;
    TYPE_PRECISION eps = EPS;

    if (LEGEND_PRINT) {
        printf("\n\n COMPLEX_3D_EX4\n");
        printf(" Forward-Backward 3D complex transform for double/single precision data\n\n");
        printf(" Configuration parameters:\n\n");
        printf(" DFTI_FORWARD_DOMAIN = DFTI_COMPLEX\n");
        printf(" DFTI_PRECISION      = DFTI_DOUBLE/DFTI_SINGLE\n");
        printf(" DFTI_DIMENSION      = 3\n");
        printf(" DFTI_LENGTHS        = {%d,%d,%d}\n", m, n, k);
        printf(" DFTI_PLACEMENT      = DFTI_NOT_INPLACE\n");
        printf(" DFTI_INPUT_STRIDES  = {0,%d,%d,%d}\n", n*k*istride, k*istride, istride);
        printf(" DFTI_OUTPUT_STRIDES = {0,%d,%d,%d}\n", n*k*ostride, k*ostride, ostride);
        printf(" DFTI_FORWARD_SCALE  = 1.0\n");
        printf(" DFTI_BACKWARD_SCALE = 1.0/(m*n*k)\n\n");
    }

    /*
    **  Allocate array for input/output and expected data
    */
    lda = m*n*k*istride*multiple;
    if (idist > m*n*k*istride) lda = idist*multiple;
    lda_out = m*n*k*ostride*multiple;
    if (odist > m*n*k*ostride) lda_out = odist*multiple;

    x_in = (fftw_complex*)fftw_malloc(2*lda*sizeof(TYPE_PRECISION));
    x_exp = (fftw_complex*)fftw_malloc(2*lda*sizeof(TYPE_PRECISION));
    x_out = (fftw_complex*)fftw_malloc(2*lda_out*sizeof(TYPE_PRECISION));

    /*
    **  Put input data and expected result
    */
#ifdef MKL_DOUBLE
    zero_init_z(x_in, lda);
    init_multiple_columns_step_z(x_in, m*n*k, multiple, idist, istride);
    cblas_zcopy(lda, x_in, 1, x_exp, 1);
#else
    zero_init_c(x_in, lda);
    init_multiple_columns_step_c(x_in, m*n*k, multiple, idist, istride);
    cblas_ccopy(lda, x_in, 1, x_exp, 1);
#endif

    /*
    **  Create FFTW plan for 3D double/single precision forward transform
    */
    my_plan = fftw3d_create_plan_specific(m, n, k, FFTW_FORWARD, FFTW_ESTIMATE,
                                          x_in, istride, x_in, ostride);
    DIE_UNLESS(my_plan);

    /*
    **  Compute DFT
    */
    fftwnd(my_plan, multiple, x_in, istride, idist, x_out, ostride, odist);

    /*
    **  Destroy FFTW plan
    */
    fftwnd_destroy_plan(my_plan);

    /*
    **  Set Scale number for Backward transform
    */
    Scale = 1.0/(TYPE_PRECISION)(m*n*k);

    /*
    **  Create FFTW plan for 3D double/single precision backward transform
    */
    my_plan = fftw3d_create_plan_specific(m, n, k, FFTW_BACKWARD, FFTW_ESTIMATE,
                                          x_out, ostride, x_in, istride);
    DIE_UNLESS(my_plan);

    /*
    **  Compute DFT
    */
    fftwnd(my_plan, multiple, x_out, ostride, odist, x_in, istride, idist);

    /*
    **  Destroy FFTW plan
    */
    fftwnd_destroy_plan(my_plan);

    /*
    **  Result scaling
    */
    scaling_multiple(x_in, Scale, m*n*k, multiple, istride, idist);

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
