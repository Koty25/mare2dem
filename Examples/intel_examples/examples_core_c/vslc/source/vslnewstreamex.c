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
!    vslNewStreamEx  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    7777777
#define N       1000

int main()
{
    VSLStreamStatePtr stream;
    VSLStreamStatePtr streamEx;
    unsigned int seed;
    unsigned int seedEx[6];
    int r[N];
    int rEx[N];
    int i, err = 0, errcode;

    /***** Initialize seeds *****/
    seed = SEED;
    seedEx[0] = SEED;
    for(i = 1; i < 6; i++){
      seedEx[i] = 1;
    }

    /***** Initialize streams *****/
    errcode = vslNewStream  ( &stream,   VSL_BRNG_MRG32K3A,    (MKL_INT)seed );
    CheckVslError( errcode );
    errcode = vslNewStreamEx( &streamEx, VSL_BRNG_MRG32K3A, 6,          seedEx );
    CheckVslError( errcode );

    /***** Call RNGs *****/
    errcode = viRngUniformBits( VSL_RNG_METHOD_UNIFORMBITS_STD, stream,   N, (unsigned int *)(r) );
    CheckVslError( errcode );
    errcode = viRngUniformBits( VSL_RNG_METHOD_UNIFORMBITS_STD, streamEx, N, (unsigned int *)(rEx) );
    CheckVslError( errcode );

    /***** Compare results *****/
    for(i = 0; i < N; i++){
      if(r[i] != rEx[i])
          err++;
    }

    /***** Printing results *****/
    printf(" Sample of vslNewStreamEx\n");
    printf(" ------------------------\n\n");
    printf(" Parameters:\n");
    printf("    seed   =   %d\n",seed);
    printf("    seedEx = { %d %d %d %d %d %d }\n\n",
      seedEx[0],seedEx[1],seedEx[2],
      seedEx[3],seedEx[4],seedEx[5]);


    printf(" Results (first 10 of 1000):\n");
    printf(" ---------------------------\n");
    for(i=0;i<10;i++) {
        printf("r[%d]=0x%08X rEx[%d]=0x%08X\n",i,r[i],i,rEx[i]);
    }

    printf("\n");
    if(err) {
        printf("Error: %d values are incorrect!\n", err);
        return 1;
    }
    else {
        printf(" Results of ordinary and extended NewStream functions are identical.\n");
    }

    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );
    errcode = vslDeleteStream( &streamEx );
    CheckVslError( errcode );

    return 0;
}
