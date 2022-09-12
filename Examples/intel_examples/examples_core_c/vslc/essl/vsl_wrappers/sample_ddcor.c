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

void ddcor(
    double h[], int inch,
    double x[], int incx,
    double y[], int incy,
    int nh, int nx, MKL_INT iy0, int ny, MKL_INT id)
{
    int status = VSL_STATUS_OK, error;
    VSLCorrTaskPtr task, task_ptr=&task;

    vsldCorrNewTask1D(task_ptr,VSL_CORR_MODE_DIRECT,nh,nx,ny);
    vslCorrSetStart(task,&iy0);
    vslCorrSetDecimation(task,&id);
    status = vsldCorrExec1D(task,h,inch,x,incx,y,incy);

    error = vslCorrDeleteTask(task_ptr);

    if (status != VSL_STATUS_OK) {
        printf("ERROR: ddcor(): bad status=%d\n",status);
        exit(1);
    }
    if (error != 0) {
        printf("ERROR: ddcor(): failed to destroy the task descriptor\n");
        exit(1);
    }
}
