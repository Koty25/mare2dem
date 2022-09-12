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
!    stream2file functions  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vsl.h"
#include "errcheck.inc"

/* Quantity of random numbers to generate */
#define N 10

static float r_orig[N], r_load[N];

int main(void)
{
    VSLStreamStatePtr stream;
    int i, errcode;

    /* Create the original stream to be saved in a file */
    errcode = vslNewStream(&stream, VSL_BRNG_R250, 7777777);
    CheckVslError( errcode );

    /* Save original stream to a file */
    errcode = vslSaveStreamF(stream, "vslstream2file.dat");
    CheckVslError( errcode );

    /* Generate random numbers using original stream */
    errcode = vsRngUniform(0, stream, N, r_orig, 0.0f, 1.0f);
    CheckVslError( errcode );

    /* Delete original stream */
    errcode = vslDeleteStream(&stream);
    CheckVslError( errcode );

    /* Load stream that is saved in a file */
    errcode = vslLoadStreamF(&stream, "vslstream2file.dat");
    CheckVslError( errcode );

    /* Generate random numbers using the stream loaded from file */
    errcode = vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, N, r_load, 0.0f, 1.0f);
    CheckVslError( errcode );

    /* Delete stream loaded from file */
    errcode = vslDeleteStream(&stream);
    CheckVslError( errcode );

    /* Compare random numbers from original and loaded stream.
       Must be identical */
    for ( i=0; i<N; i++ )
    {
        printf("r_orig[%d]=%f\tr_load[%d]=%f\n", i, r_orig[i], i, r_load[i]);
        if ( r_orig[i] != r_load[i] )
        {
            /* Here if results are not identical */
            printf("Error: Loaded stream differs from original stream.\n");
            return 1;
        }
    }

    /* Here if results are identical */
    printf("PASS: Loaded stream identical with original stream.\n");

    return 0;
}
