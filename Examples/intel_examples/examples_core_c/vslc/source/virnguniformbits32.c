/*******************************************************************************
* Copyright 2003-2020 Intel Corporation.
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
!    viRngUniformBits32  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    1
#define BRNG    VSL_BRNG_MCG59
#define METHOD  VSL_RNG_METHOD_UNIFORMBITS32_STD
#define N       1000
#define NN      10

int main()
{
    int r[N];
    VSLStreamStatePtr stream;
    int i, errcode;

    /***** Initialize *****/
    errcode = vslNewStream( &stream, BRNG,  SEED );
    CheckVslError( errcode );

    /***** Call RNG *****/
    errcode = viRngUniformBits32( METHOD, stream, N, (unsigned int *)(r) );
    CheckVslError( errcode );

    /***** Printing results *****/
    printf("Sample of viRngUniformBits32.\n");
    printf("---------------------------\n\n");

    printf("Results (first 10 of 1000):\n");
    printf("---------------------------\n");
    for(i=0;i<NN;i++) {
        printf("r[%d]=%d\n",i,r[i]);
    }

    printf("\n");

    /***** Deinitialize *****/
    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    return 0;
}
