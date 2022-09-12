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
!    stream2memory functions  Example Program Text
!******************************************************************************/

#include <stdio.h>
#include "mkl_vsl.h"
#include "errcheck.inc"

/* Quantity of random numbers to generate */
#define DN 101

#define N   10

#define SEED       7777777
#define BRNG       VSL_BRNG_SFMT19937
#define METHOD     VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2

int main(void)
{
    double rd_orig[DN], rd_load[DN];
    int size;
    char *membuf;

    VSLStreamStatePtr stream_orig, stream_load;
    int i, errcode;

    /* Create the original stream to be saved in a file */
    errcode = vslNewStream( &stream_orig, BRNG , SEED );
    CheckVslError( errcode );

    /* Generate random numbers using original stream */
    errcode = vdRngGaussian( METHOD, stream_orig, DN, rd_orig, 0.0, 1.0);
    CheckVslError( errcode );

    /* Generation of Gaussian random numbers */
    /* Compute memory size necessary for random stream descriptive data */
    size = vslGetStreamSize( stream_orig );
    membuf = (char*)malloc( sizeof(char) * size );
    if ( membuf == 0 )
    {
        printf("memeory allocation error\n");
        return 1;
    }
    errcode = vslSaveStreamM( stream_orig, membuf );
    CheckVslError( errcode );

    /* Generate random numbers using original stream */
    errcode = vdRngGaussian( METHOD, stream_orig, DN, rd_orig, 0.0, 1.0);
    CheckVslError( errcode );

    /* Load stream from the memeory */
    errcode = vslLoadStreamM( &stream_load, membuf );
    CheckVslError( errcode );

    /* Generate random numbers using loaded stream */
    errcode = vdRngGaussian( METHOD, stream_load, DN, rd_load, 0.0, 1.0);
    CheckVslError( errcode );

    /* Compare random numbers from original and loaded stream.
       Must be identical */
    printf("Gaussian numbers:\n");
    for ( i=0; i<N; i++ )
    {
        printf("rd_orig[%d]=%f\trd_load[%d]=%f\n", i, rd_orig[i], i, rd_load[i]);
    }

    for ( i=0; i<DN; i++ )
    {
        if ( rd_orig[i] != rd_load[i] )
        {
            /* Here if results are not identical */
            printf("Error: Loaded stream differs from original stream.\n");
            return 1;
        }
    }

    free( membuf );

    /* Delete stream loaded from memory */
    errcode = vslDeleteStream(&stream_load);
    CheckVslError( errcode );

    /* Delete original stream */
    errcode = vslDeleteStream(&stream_orig);
    CheckVslError( errcode );

    /* Here if results are identical */
    printf("PASS: Loaded stream identical with original stream.\n");

    return 0;
}
