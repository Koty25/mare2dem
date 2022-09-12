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
!  Example of 1-dimension linear convolution operation on double precision data.
!*******************************************************************************/

#include <math.h>
#include <stdio.h>

#include "mkl_vsl.h"

int main()
{
    VSLConvTaskPtr task;
    double x[4]={1,2,3,4};
    double y[8]={11,12,13,14,15,16,17,18};
    double z[11]={0,0,0,0,0,0,0,0,0,0,0};
    double e[11]={11,34,70,120,130,140,150,160,151,122,72};
    MKL_INT xshape=4, yshape=8, zshape=11;
    int status, i;

    int mode = VSL_CONV_MODE_AUTO;

    /*
    *  Create task descriptor (create descriptor of problem)
    */
    status = vsldConvNewTask1D(&task,mode,xshape,yshape,zshape);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: creation of job failed, exit with %d\n", status);
        return 1;
    }

    /*
    *  Execute task (Calculate 1 dimension convolution of two arrays)
    */
    status = vsldConvExec1D(task,x,1,y,1,z,1);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: job status bad, exit with %d\n", status);
        return 1;
    }

    /*
    *  Delete task object (delete descriptor of problem)
    */
    status = vslConvDeleteTask(&task);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: failed to delete task object, exit with %d\n", status);
        return 1;
    }

    /*
    * Check resulst for correctness:
    */
    for (i=0; i<zshape; i++)
        if (fabs(z[i]-e[i]) > fabs(e[i])*1e-10) {
            printf("ERROR: wrong results:\n");
            printf("    z[%2d]: %lg\n",i,z[i]);
            printf(" expected: %lg\n",e[i]);
            printf("EXAMPLE FAILED\n");
            return 1;
        }

    printf("EXAMPLE PASSED\n");
    return 0;
}
