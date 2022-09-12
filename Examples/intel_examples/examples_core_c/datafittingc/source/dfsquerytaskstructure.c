/*******************************************************************************
* Copyright 2010-2020 Intel Corporation.
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
!   Querying Data Fitting task parameters Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N                     10    // number of breakpoints
#define NY                     2    // number of datasets to interpolate

int main()
{
    DFTaskPtr task;                 // Data Fitting task descriptor
    MKL_INT nx;                     // number of break points
    MKL_INT xhint;                  // additional info about break points
    MKL_INT ny;                     // number of functions
    MKL_INT yhint;                  // functions storage format

    float x[] = {-1.0f, 1.0f};      // break points
    float y[N*NY];                  // function values

    MKL_INT nx_ret;                 // parameters to be queried
    float *y_ret;
    float *y1_ret;

    int errnums, errcode = 0;

    /***** Initializing parameters for Data Fitting task *****/

    /***** Parameters describing interpolation interval *****/
    nx        = N;
    xhint     = DF_UNIFORM_PARTITION;

    /***** Parameters describing functions *****/
    ny         = NY;
    yhint      = DF_MATRIX_STORAGE_ROWS;

    /***** Create Data Fitting task *****/
    errcode = dfsNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Query Data Fitting task parameters *****/

    /***** Query value *****/
    errcode = dfiQueryVal( task, DF_NX, &nx_ret );
    CheckDfError(errcode);

    /***** Query pointer *****/
    errcode = dfsQueryPtr( task, DF_Y, &y_ret );
    CheckDfError(errcode);

    /***** Query pointer by index *****/
    errcode = dfsQueryIdxPtr( task, DF_Y, 1, &y1_ret );
    CheckDfError(errcode);

    /***** Check requested parameters *****/
    errnums = 0;
    if ( nx    != nx_ret ) errnums++;
    if ( y     != y_ret )  errnums++;
    if ( y + N != y1_ret ) errnums++;

    /***** Print results *****/
    printf("                                          Expected");
    printf("              Obtained\n");
    printf("Number of break points          : %16d      %16d\n",
        (int)nx, (int)nx_ret);
    printf("Address of function Y           : %16p      %16p\n", y, y_ret);
    printf("Address of 1-th coordinate of Y : %16p      %16p\n", y+N, y1_ret);

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Print summary of the test *****/
    if (errnums != 0)
    {
        printf("\n\nError: Not all requested parameters are correct\n");
        return 1;
    }
    else
    {
        printf("\n\nAll requested parameters are correct\n");
    }

    return 0;
}
