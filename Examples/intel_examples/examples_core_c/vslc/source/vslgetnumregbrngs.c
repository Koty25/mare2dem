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
!    vslGetNumRegBrngs  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl_vsl.h"
#include "errcheck.inc"

#define BRNGS   16

int main()
{
    int brngsExp = BRNGS;
    int brngsObt = 0;

    /***** Get number of registered BRNGs *****/
    brngsObt = vslGetNumRegBrngs ();

    /***** Printing results *****/
    printf(" Sample of vslGetNumRegBrngs\n");
    printf(" -------------------------------\n\n");

    if(brngsObt != brngsExp) {
        printf(" Error: returned value %d is incorrect (expected %d)!\n", brngsObt,brngsExp);
        return 1;
    }
    else {
        printf(" Returned %d as expected\n",brngsObt);
    }

    return 0;
}
