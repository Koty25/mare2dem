/*******************************************************************************
* Copyright 2019-2020 Intel Corporation.
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

#ifndef _OPENMP_COMMON_H_
#define _OPENMP_COMMON_H_

#define MAX(a,b) ((a) > (b)) ? (a) : (b)

static inline double rand_double_scalar() {
    return ((double) rand () / (double) RAND_MAX) - 0.5; 
}

static inline float rand_single_scalar() {
    return ((float) rand () / (float) RAND_MAX) - 0.5; 
}

static inline MKL_Complex16 rand_double_complex_scalar() {
    MKL_Complex16 res;
    res.real = ((double) rand () / (double) RAND_MAX) - 0.5; 
    res.imag = ((double) rand () / (double) RAND_MAX) - 0.5; ;
    return res;
}

static inline MKL_Complex8 rand_single_complex_scalar() {
    MKL_Complex8 res;
    res.real = ((float) rand () / (float) RAND_MAX) - 0.5; 
    res.imag = ((float) rand () / (float) RAND_MAX) - 0.5; ;
    return res;
}

static inline void free_double_matrices(double **array, MKL_INT size) {
    for (MKL_INT i = 0; i < size; i++) {
        mkl_free(array[i]);
    }
}

static inline void free_single_matrices(float **array, MKL_INT size) {
    for (MKL_INT i = 0; i < size; i++) {
        mkl_free(array[i]);
    }
}

static inline void free_double_complex_matrices(MKL_Complex16 **array, MKL_INT size) {
    for (MKL_INT i = 0; i < size; i++) {
        mkl_free(array[i]);
    }
}

static inline void free_single_complex_matrices(MKL_Complex8 **array, MKL_INT size) {
    for (MKL_INT i = 0; i < size; i++) {
        mkl_free(array[i]);
    }
}

static inline MKL_INT rand_int(MKL_INT min, MKL_INT max) {
    MKL_INT res = min + (rand() % (max - min + 1));
    return res;
}

static inline void init_double_array(MKL_INT n, double *array, MKL_INT do_rand) {
    MKL_INT i;
    for (i = 0; i < n; i++) {
        if (do_rand) {
            array[i] = rand() / (double) RAND_MAX - .5;
        }
        else {
            array[i] = (double) (i + 1);
        }
    }
}

static inline void init_single_array(MKL_INT n, float *array, MKL_INT do_rand) {
    MKL_INT i;
    for (i = 0; i < n; i++) {
        if (do_rand) {
            array[i] = rand() / (float) RAND_MAX - .5;
        }
        else {
            array[i] = (float) (i + 1);
        }
    }
}

static inline void init_double_complex_array(MKL_INT n, MKL_Complex16 *array, MKL_INT do_rand) {
    MKL_INT i;
    for (i = 0; i < n; i++) {
        if (do_rand) {
            array[i].real = rand() / (double) RAND_MAX - .5;
            array[i].imag = rand() / (double) RAND_MAX - .5;
        }
        else {
            array[i].real = (double) (i + 1);
            array[i].imag = (double) (i + 1);
        }
    }
}

static inline void init_single_complex_array(MKL_INT n, MKL_Complex8 *array, MKL_INT do_rand) {
    MKL_INT i;
    for (i = 0; i < n; i++) {
        if (do_rand) {
            array[i].real = rand() / (float) RAND_MAX - .5;
            array[i].imag = rand() / (float) RAND_MAX - .5;
        }
        else {
            array[i].real = (float) (i + 1);
            array[i].imag = (float) (i + 1);
        }
    }
}

#endif
