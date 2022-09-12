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
!    vslGetStreamStateBrng  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define SEED    7777777

int main()
{
    VSLStreamStatePtr stream;
    unsigned int seed;
    int errcode;
    MKL_INT brngExp = VSL_BRNG_WH+127;
    int brngObt = 0;

    /***** Initialize seed *****/
    seed = SEED;

    /***** Initialize streams *****/
    errcode = vslNewStream  ( &stream, brngExp, (MKL_INT)seed );
    CheckVslError( errcode );

    /***** Get BRNG number *****/
    brngObt = vslGetStreamStateBrng ( stream );

    /***** Printing results *****/
    printf(" Sample of vslGetStreamStateBrng\n");
    printf(" -------------------------------\n\n");
    printf(" Parameters:\n");
    printf("    seed = %d\n",seed);
    printf("    brng = %d\n\n", (int)brngExp);

    if(brngObt != brngExp) {
        printf(" Error: returned value %d is incorrect (expected %d)!\n", brngObt, (int)brngExp);
        return 1;
    }
    else {
        printf(" Returned %d as expected\n",brngObt);
    }

    errcode = vslDeleteStream( &stream );
    CheckVslError( errcode );

    return 0;
}
