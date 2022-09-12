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
!    vslConvDeleteTask  Example Program Text
!******************************************************************************/

#include "mkl_vsl.h"

#include <stdio.h>

int main()
{
    VSLConvTaskPtr task;
    int ok,error;
    MKL_INT mode,rank,xshape,yshape,zshape;

    ok = 1;
    printf("EXAMPLE deleting a task\n");

    mode = VSL_CONV_MODE_AUTO;
    rank = 1;
    xshape = 100;
    yshape = 1000;
    zshape = (xshape-1) + (yshape-1) + 1;
    vslsConvNewTask(&task,mode,rank,&xshape,&yshape,&zshape);

    error = vslConvDeleteTask(&task);

    if (error != 0) {
        printf("ERROR: error code: %d\n",error);
        ok = 0;
    }
    if (task != NULL) {
        printf("ERROR: task pointer isn't cleaned!\n");
        ok = 0;
    }

    printf("EXAMPLE %s\n", ok? "PASSED": "FAILED");
    return !ok;
}
