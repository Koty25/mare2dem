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

void scond(
    float h[], int inch,
    float x[], int incx,
    float y[], int incy,
    int nh, int nx, MKL_INT iy0, int ny)
{
    int status = VSL_STATUS_OK, error;
    VSLConvTaskPtr task, task_ptr=&task;
    MKL_INT *start = &iy0; /* simulate start[] array */

    vslsConvNewTask1D(task_ptr,VSL_CONV_MODE_DIRECT,nh,nx,ny);
    vslConvSetStart(task,start);
    status = vslsConvExec1D(task,h,inch,x,incx,y,incy);

    error = vslConvDeleteTask(task_ptr);

    if (status != VSL_STATUS_OK) {
        printf("ERROR: scond(): bad status=%d\n",status);
        exit(1);
    }
    if (error != 0) {
        printf("ERROR: scond(): failed to destroy the task descriptor\n");
        exit(1);
    }
}
