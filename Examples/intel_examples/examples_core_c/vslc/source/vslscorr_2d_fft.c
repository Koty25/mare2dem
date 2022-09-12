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
!  Example of 2-dimension linear correlation operation on single precision data.
!*******************************************************************************/

#include <math.h>
#include <stdio.h>

#include "mkl_vsl.h"

int main()
{
    VSLCorrTaskPtr task;
    float x[3*2]={1,1,1, 1,1,1};
    float y[4*3]={1,1,1,1, 1,1,1,1, 1,1,1,1};
    float z[6*4]={0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0, 0,0,0,0,0,0};
    float e[6*4]={1,2,3,3,2,1, 2,4,6,6,4,2, 2,4,6,6,4,2, 1,2,3,3,2,1};
    MKL_INT xshape[2]={3,2}, yshape[2]={4,3}, zshape[2]={6,4};
    MKL_INT rank=2;
    int status,i,j;

    int mode=VSL_CORR_MODE_FFT;

    /*
    *  Create task descriptor (create descriptor of problem)
    */
    status = vslsCorrNewTask(&task,mode,rank,xshape,yshape,zshape);
    if( status != VSL_STATUS_OK ){
        printf("ERROR: creation of job failed, exit with %d\n", status);
        return 1;
    }

    /*
    *  Execute task (Calculate 2 dimension correlation of two arrays)
    */
    status = vslsCorrExec(task,x,NULL,y,NULL,z,NULL);
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

    for (j=0; j<zshape[1]; j++)
        for (i=0; i<zshape[0]; i++) {
            float zij = z[i + zshape[0]*j];
            float eij = e[i + zshape[0]*j];
            if (fabs(zij-eij) > fabs(eij)*1e-5) {
                printf("ERROR: wrong results:\n");
                printf("    z[%2d,%2d]: %g\n",i,j,zij);
                printf("    expected: %g\n",eij);
                printf("EXAMPLE FAILED\n");
                return 1;
            }
        }

    printf("EXAMPLE PASSED\n");
    return 0;
}
