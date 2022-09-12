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

#include "mkl_vsl.h"

#include <stdio.h>
#include <stdlib.h>

void scorf(
    int init,
    float h[], int inc1h,
    float x[], int inc1x, int inc2x,
    float y[], int inc1y, int inc2y,
    int nh, int nx, int m, MKL_INT iy0, int ny,
    void* aux1, int naux1, void* aux2, int naux2)
{
    int status = VSL_STATUS_OK ,error,i;
    VSLCorrTaskPtr task;

    if (init != 0)
        return; /* ignore aux1, aux2 */

    vslsCorrNewTaskX1D(&task,
        VSL_CORR_MODE_FFT,nh,nx,ny,
        h,inc1h);
    vslCorrSetStart(task, &iy0);

    /* task is implicitly committed at i==0 */
    for (i=0; i<m; i++) {
        float* xi = &x[inc2x * i];
        float* yi = &y[inc2y * i];
        status = vslsCorrExecX1D(task,
            xi,inc1x, yi,inc1y);
        /* check status later */
    }

    error = vslCorrDeleteTask(&task);

    if (status != VSL_STATUS_OK) {
        printf("ERROR: scorf(): bad status=%d\n",status);
        exit(1);
    }
    if (error != 0) {
        printf("ERROR: scorf(): failed to destroy the task descriptor\n");
        exit(1);
    }
}
