/*******************************************************************************
* Copyright 2001-2020 Intel Corporation.
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
!  Example of 1-dimension linear correlation operation on single precision data.
!*******************************************************************************/

#include <math.h>
#include <stdio.h>

#include "mkl_vsl.h"

int main()
{
    VSLCorrTaskPtr task;
    float x[4]={1,2,3,4};
    float y[8]={11,12,13,14,15,16,17,18};
    float z[11]={0,0,0,0,0,0,0,0,0,0,0};
    float e[11]={44,81,110,130,140,150,160,170,104,53,18};
    MKL_INT xshape=4, yshape=8, zshape=11;
    int status,i;

    int mode = VSL_CORR_MODE_AUTO;

    /*
    *  Create task descriptor (create descriptor of problem)
    */
    status = vslsCorrNewTask1D(&task,mode,xshape,yshape,zshape);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: creation of job failed, exit with %d\n", status);
        return 1;
    }

    /*
    *  Execute task (Calculate correlation of two arrays)
    */
    status = vslsCorrExec1D(task,x,1,y,1,z,1);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: job status bad, exit with %d\n", status);
        return 1;
    }

    /*
    *  Delete task object (delete descriptor of problem)
    */
    status = vslCorrDeleteTask(&task);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: failed to delete task object, exit with %d\n", status);
        return 1;
    }

    /*
    * Check resulst for correctness:
    */
    for (i=0; i<zshape; i++)
        if (fabs(z[i]-e[i]) > fabs(e[i])*1e-5) {
            printf("ERROR: wrong results:\n");
            printf("    z[%2d]: %g\n",i,z[i]);
            printf(" expected: %g\n",e[i]);
            printf("EXAMPLE FAILED\n");
            return 1;
        }

    printf("EXAMPLE PASSED\n");
    return 0;
}
