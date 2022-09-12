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
!    vsldCorrExecX1D  Example Program Text
!******************************************************************************/

#include "mkl_vsl.h"

#include <stdio.h>

#define XSHAPE 100
#define YSHAPE 1000
#define ZSHAPE (XSHAPE-1)+(YSHAPE-1)+1

int main()
{
    VSLCorrTaskPtr task;
    MKL_INT mode,xshape,yshape,zshape;
    static double x[XSHAPE],y[YSHAPE],z[ZSHAPE];
    MKL_INT xstride=1,ystride=1,zstride=1;
    int status,ok,i;

    xshape=XSHAPE;
    yshape=YSHAPE;
    zshape=ZSHAPE;

    for (i=0; i<xshape; i++)
        x[i] = 0;
    for (i=0; i<yshape; i++)
        y[i] = 0;

    ok = 1;
    printf("EXAMPLE executing a correlation task\n");

    mode = VSL_CORR_MODE_AUTO;
    vsldCorrNewTaskX1D(&task,mode,xshape,yshape,zshape,
        x,xstride);

    status = vsldCorrExecX1D(task,y,ystride,z,zstride);

    if (status != VSL_STATUS_OK) {
        printf("ERROR: bad status: %d\n",status);
        ok = 0;
    }

    for (i=0; i<zshape; i++)
        if (z[i] != 0) {
            printf("ERROR: wrong result: z[%d]=%lg\n",i,z[i]);
            ok = 0;
        }

    printf("EXAMPLE %s\n", ok? "PASSED": "FAILED");
    return !ok;
}
