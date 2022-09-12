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
!    vsRngGaussian  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include <math.h>

#include "mkl.h"
#include "mkl_omp_offload.h"

#define SEED    1
#define BRNG    VSL_BRNG_MRG32K3A
#define METHOD  VSL_RNG_METHOD_GAUSSIAN_ICDF
#define N       1000
#define NN      10

int main() {
    int dnum = 0;

    float* r = (float*)mkl_malloc((N) * sizeof(float), 64);
    if(r == NULL) {
        printf("Cannot allocate memory\n");
        return 1;
    }

    VSLStreamStatePtr stream;
    int i;
    float mean = 0.0f, stddev = 1.0f;

    // initialize Basic Random Number Generator
    vslNewStream(&stream, BRNG,  SEED);

#pragma omp target data map(tofrom:r[0:N]) device(dnum)
    {
// run RNG on gpu, use standard oneMKL interface within a variant dispatch construct
#pragma omp target variant dispatch device(dnum) use_device_ptr(r)
        {
            vsRngGaussian(METHOD, stream, N, r, mean, stddev);
        }
    }

    double tM,tD,tQ,tD2;
    double sM,sD;
    double sum, sum2;
    double n,s;
    double DeltaM,DeltaD;

    // theoretical moments
    tM = mean;
    tD = stddev * stddev;
    tQ = 720.0 * stddev * stddev * stddev * stddev;

    // sample moments
    sum = 0.0;
    sum2 = 0.0;
    for(i = 0; i < N; i++) {
        sum += (double)r[i];
        sum2 += (double)r[i] * (double)r[i];
    }
    sM = sum / ((double)N);
    sD = sum2 / (double)N - (sM * sM);

    // comparison of theoretical and sample moments
    n = (double)N;
    tD2 = tD * tD;
    s = ((tQ-tD2) / n) - (2 * (tQ - 2 * tD2) / (n * n)) + ((tQ - 3 * tD2) / (n * n * n));

    DeltaM = (tM - sM) / sqrt(tD / n);
    DeltaD = (tD - sD) / sqrt(s);

    // printing results
    printf("Sample of vsRngGaussian.\n");
    printf("------------------------\n\n");
    printf("Parameters:\n");
    printf("    mean=%.4f\n", mean);
    printf("    stddev=%.4f\n\n", stddev);

    printf("Results (first 10 of 1000):\n");
    printf("---------------------------\n");
    for(i = 0; i < NN; i++) {
        printf("r[%d]=%.4f\n", i, r[i]);
    }

    printf("\n");
    if(fabs(DeltaM) > 3.0 || fabs(DeltaD) > 3.0) {
        printf("Error: sample moments (mean=%.2f, variance=%.2f) disagree with theory (mean=%.2f, variance=%.2f).\n", sM, sD, tM, tD);
        return 1;
    }
    else {
        printf("Sample moments (mean=%.2f, variance=%.2f) agree with theory (mean=%.2f, variance=%.2f).\n", sM, sD, tM, tD);
    }

    mkl_free(r);

    // deinitialize
    vslDeleteStream(&stream);

    return 0;
}
