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
!    vslConvSetStart  Example Program Text
!******************************************************************************/

#include "mkl_vsl.h"

#include <stdio.h>

int main()
{
    VSLConvTaskPtr task;
    int status,ok;
    MKL_INT mode,rank,xshape,yshape,zshape,start;

    ok = 1;
    printf("EXAMPLE setting start parameter\n");

    mode = VSL_CONV_MODE_AUTO;
    rank = 1;
    xshape = 100;
    yshape = 1000;
    zshape = (xshape-1) + (yshape-1) + 1;
    vslsConvNewTask(&task,mode,rank,&xshape,&yshape,&zshape);

    start = 0;
    status = vslConvSetStart(task,&start);

    if (status != VSL_STATUS_OK) {
        printf("ERROR: bad status: %d\n",status);
        ok = 0;
    }

    printf("EXAMPLE %s\n", ok? "PASSED": "FAILED");
    return !ok;
}
