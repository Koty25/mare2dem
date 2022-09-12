/*******************************************************************************
* Copyright 2020 Intel Corporation.
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
!    vsRngUniform  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "mkl_omp_offload.h"

#define SEED    1
#define BRNG    VSL_BRNG_PHILOX4X32X10
#define METHOD  VSL_RNG_METHOD_UNIFORM_STD
#define N       1000
#define NN      10

int main() {
    int dnum = 0;

    float* r_host = (float*)mkl_malloc((N) * sizeof(float), 64);
    float* r_dev = (float*)mkl_malloc((N) * sizeof(float), 64);
    if((r_host == NULL) || (r_dev == NULL)) {
        printf("Cannot allocate memory\n");
        return 1;
    }

    VSLStreamStatePtr stream_host;
    VSLStreamStatePtr stream_dev;
    int i;
    float a = 0.0f, b = 1.0f;

    // initialize Basic Random Number Generators
    vslNewStream(&stream_host, BRNG,  SEED);
    vslNewStream(&stream_dev, BRNG,  SEED);

    // call RNG on host
    vsRngUniform(METHOD, stream_host, N, r_host, a, b);

#pragma omp target data map(tofrom:r_dev[0:N]) device(dnum)
    {
// run RNG on gpu, use standard oneMKL interface within a variant dispatch construct
#pragma omp target variant dispatch device(dnum) use_device_ptr(r_dev)
        {
            vsRngUniform(METHOD, stream_dev, N, r_dev, a, b);
        }
    }

    // comparison of host and device results
    int err = 0;
    for(i = 0; i < N; i++) {
        if(r_dev[i] != r_host[i]) {
            printf("Error in %d element r_dev = %.4f r_host = %.4f\n", i, r_dev[i], r_host[i]);
            err++;
        }
    }

    // printing results
    printf("Sample of vsRngUniform.\n");
    printf("-----------------------\n\n");
    printf("Parameters:\n");
    printf("    a=%.4f\n",a);
    printf("    b=%.4f\n\n",b);

    printf("Results (first 10 of 1000):\n");
    printf("---------------------------\n");
    for(i = 0; i < NN; i++) {
        printf("r_dev[%d]=%.4f\n", i, r_dev[i]);
    }

    mkl_free(r_host);
    mkl_free(r_dev);

    // deinitialize
    vslDeleteStream(&stream_host);
    vslDeleteStream(&stream_dev);

    return err;
}
