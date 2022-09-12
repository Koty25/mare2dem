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
!    vslCopyStreamState  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    7777777
#define N       1000

int main()
{
    VSLStreamStatePtr stream;
    VSLStreamStatePtr streamCpy;
    unsigned int seedCpy;
    unsigned int seed[6];
    int r[N];
    int rCpy[N];
    int i, err = 0, errcode;

    /***** Initialize seeds *****/
    seedCpy = 1;
    seed[0] = SEED;
    for(i = 1; i < 6; i++){
      seed[i] = seed[i-1]+11;
    }

    /***** Initialize streams *****/
    errcode = vslNewStreamEx( &stream,    VSL_BRNG_MRG32K3A, 6,          seed );
    CheckVslError( errcode );
    errcode = vslNewStream  ( &streamCpy, VSL_BRNG_MRG32K3A,    (MKL_INT)seedCpy );
    CheckVslError( errcode );
    errcode = vslCopyStreamState ( streamCpy, stream );
    CheckVslError( errcode );


    /***** Call RNGs *****/
    errcode = viRngUniformBits( VSL_RNG_METHOD_UNIFORMBITS_STD, stream,   N, (unsigned int *)(r) );
    CheckVslError( errcode );
    errcode = viRngUniformBits( VSL_RNG_METHOD_UNIFORMBITS_STD, streamCpy, N, (unsigned int *)(rCpy) );
    CheckVslError( errcode );

    /***** Compare results *****/
    for(i = 0; i < N; i++){
      if(r[i] != rCpy[i])
          err++;
    }

    /***** Printing results *****/
    printf(" Sample of vslCopyStreamState\n");
    printf(" ----------------------------\n\n");
    printf(" Parameters:\n");
    printf("    seed    = { %d %d %d %d %d %d }\n",
      seed[0],seed[1],seed[2],
      seed[3],seed[4],seed[5]);
    printf("    seedCpy =   %d\n\n",seedCpy);


    printf(" Results (first 10 of 1000):\n");
    printf(" ---------------------------\n");
    for(i=0;i<10;i++) {
        printf("r[%d]=0x%08X rCpy[%d]=0x%08X\n",i,r[i],i,rCpy[i]);
    }

    printf("\n");
    if(err) {
        printf("Error: %d values are incorrect!\n", err);
        return 1;
    }
    else {
        printf(" Results of original stream and its copy are identical.\n");
    }

    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );
    errcode = vslDeleteStream( &streamCpy );
    CheckVslError( errcode );

    return 0;
}
