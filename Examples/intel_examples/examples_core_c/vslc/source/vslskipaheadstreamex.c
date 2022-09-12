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

/*
!  Content:
!    vslSkipaheadStreamEx  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl_vsl.h"
#include "errcheck.inc"
#include "math.h"

#define SEED    7777777
#define N       1000
#define NPARAMS 2

int main()
{
    VSLStreamStatePtr stream;
    VSLStreamStatePtr streamS;

    unsigned int r[N];

    unsigned int rS[N];

    int seed = SEED, i=0, j, err = 0, errcode;

    /*  To skip 2^76 elements in the random stream SkipAheadStream(nskip) function should be called 2^14 times
        with nskip equal to 2^62*/
    MKL_INT64 nskip = pow(2,62);
    int skipTimes = pow(2,14);

    /*  To skip 2^76 elements in the random stream SkipAheadStreamEx(nskip) function should be called
        with nskip represented as
            nskip = 2^76 = 0 + 2^12 * 2^64
        In general case
            nskip = params[0] + params[1] * 2^64 + params[2] * 2^128 + ... */
    MKL_UINT64 params[NPARAMS];
    params[0] = 0;
    params[1] = pow(2,12);

    /****** Create main stream *********/
    errcode = vslNewStream  ( &stream,   VSL_BRNG_MRG32K3A,  (MKL_INT)seed );
    CheckVslError( errcode );

    /* Create streamS as copy of the main one */
    errcode = vslCopyStream( &streamS, stream );
    CheckVslError( errcode );

    /* Apply vslSkipAheadStreamEx to the streamS with the nskip = 0 + 2^12 * 2^64 = 2^76 */
    errcode = vslSkipAheadStreamEx( streamS, NPARAMS, params );
    CheckVslError( errcode );

    /* Apply vslSkipAheadStream to the stream 2^14 times with the nskip = 2^62 */
    for (j=0; j<skipTimes; j++){
        errcode = vslSkipAheadStream( stream, nskip);
        CheckVslError( errcode );
    }

    /**** Generate random numbers for skipahead stream  ****/
    errcode = viRngUniformBits( VSL_RNG_METHOD_UNIFORMBITS_STD, stream, N, r );
    CheckVslError( errcode );

    /* Generate random numbers for skipaheadex stream  */
    errcode = viRngUniformBits( VSL_RNG_METHOD_UNIFORMBITS_STD, streamS, N, rS );
    CheckVslError( errcode );

    /***** Compare results *****/
    for ( i=0; i<N; i++ ){
      if(r[i] != rS[i])
          err++;
    }

    /***** Printing results *****/
    printf(" Sample of vslSkipaheadStreamEx\n");
    printf(" ----------------------------\n\n");
    printf(" Parameters:\n");
    printf("    seed   =   %d\n\n", seed);


    printf(" Results (first 10 of 1000):\n");
    printf(" ---------------------------\n");
    for(i=0;i<10;i++) {
        printf("r[%d]=0x%08X rS[%d]=0x%08X\n",i,r[i],i,rS[i]);
    }

    printf("\n");
    if(err) {
        printf(" Error: %d values are incorrect!\n", err);
        return 1;
    }
    else {
        printf(" Results of SkipAhead and SkipaheadEx streams are identical.\n");
    }

    errcode = vslDeleteStream( &stream );
    errcode = vslDeleteStream( &streamS );

    return 0;
}
