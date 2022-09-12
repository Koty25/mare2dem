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
!    vslLeapfrogStream  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    7777777
#define N       1000
#define S       10
#define NS      100


int main()
{
    VSLStreamStatePtr stream;
    VSLStreamStatePtr streamL[S];
    int r [N];
    int rL[N];
    int seed = SEED, i, j, k, err = 0, errcode;

    /****** Create main stream *********/
    errcode = vslNewStream  ( &stream,   VSL_BRNG_MCG31,  (MKL_INT)seed );
    CheckVslError( errcode );
    /* Create leapfrog streams as copies of the main one */
    for(i=0;i<S;i++)
    {
        errcode = vslCopyStream( &streamL[i], stream );
        CheckVslError( errcode );
        errcode = vslLeapfrogStream( streamL[i], (MKL_INT)i, S );
        CheckVslError( errcode );
    }

    /**** Generate random numbers for main stream  ****/
    errcode = viRngUniformBits( 0, stream, N, (unsigned int *)(r) );
    CheckVslError( errcode );
    /* Generate random numbers for leapfrog streams  */
    for(i=0;i<S;i++)
    {
        errcode = viRngUniformBits( VSL_RNG_METHOD_UNIFORMBITS_STD, streamL[i], NS, (unsigned int *)(&(rL[i*NS])) );
        CheckVslError( errcode );
    }

    /***** Compare results *****/
    for ( j=0, i=0; i<NS; i++ )
    {
       for ( k=0; k<S; k++ )
       {
      if(r[j++] != rL[k*NS+i])
          err++;
       }
    }

    /***** Printing results *****/
    printf(" Sample of vslLeapfrogStream\n");
    printf(" ---------------------------\n\n");
    printf(" Parameters:\n");
    printf("    seed   =   %d\n\n",seed);


    printf(" Results (first 10 of 1000):\n");
    printf(" ---------------------------\n");
    for(i=0;i<10;i++) {
        printf("r[%d]=0x%08X rL[%3d]=0x%08X\n",i,r[i],i*NS,rL[i*NS]);
    }

    printf("\n");
    if(err) {
        printf("Error: %d values are incorrect!\n", err);
        return 1;
    }
    else {
        printf(" Results of ordinary and Leapfrog streams are identical.\n");
    }

    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );
    for(i=0;i<S;i++)
    {
        errcode = vslDeleteStream( &streamL[i] );
        CheckVslError( errcode );
    }

    return 0;
}
