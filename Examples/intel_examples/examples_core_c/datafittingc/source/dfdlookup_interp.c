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
!    Interpolation with Look up method
!    Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N              7 // number of break points
#define NY             1 // number of functions
#define NSITE          N // number of sites

#define NDORDER        1 // size of array describing derivative orders
                         // to compute
#define NDER     NDORDER // number of derivatives to compute

#define LEFT_LIMIT   1.0 // left limit of interpolation interval
#define RIGHT_LIMIT  3.0 // right limit of interpolation interval
#define FREQ         0.5


int main()
{
    DFTaskPtr task;                     // Data Fitting task descriptor
    MKL_INT nx;                         // number of break points
    MKL_INT xhint;                      // additional info about break points
    MKL_INT ny;                         // number of functions
    MKL_INT yhint;                      // additional info about function

    double x[N];                        // array of break points
    double y[NY*N];                     // function values

    MKL_INT nsite;                      // total number of interpolation sites
    MKL_INT sitehint;                   // additional info about interpolation
    double site[NSITE];                 // array of interpolation sites

    MKL_INT ndorder;                    // size of array describing derivative
                                        // orders
    MKL_INT dorder[] = { 1 };           // only value to calculate
                                        // will be computed

    MKL_INT rhint;                      // interpolation results storage format
    double r[NSITE];                    // spline evaluation results

    MKL_INT stype, sorder;

    double *datahint;                   // additional info about structure
                                        // of arrays x and y

    double left = LEFT_LIMIT, right = RIGHT_LIMIT;
    double freq =FREQ;

    int i, j, errcode = 0;
    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/
    sorder = DF_PP_STD;
    stype  = DF_LOOKUP_INTERPOLANT;

    /***** Additional info about structure of arrays x and y *****/
    /* No additional info is provided */
    datahint = 0;

    /***** Parameters describing interpolation interval *****/
    nx          = N;
    xhint       = DF_NON_UNIFORM_PARTITION;

    /***** Parameters describing function *****/
    ny          = NY;
    yhint       = DF_NO_HINT;

    /***** Additional info about structure of arrays x and y *****/
    /* No additional info is provided */
    datahint = 0;

    /***** Parameters describing interpolation sites *****/
    nsite      = NSITE;
    sitehint   = DF_NON_UNIFORM_PARTITION;

    /**** Parameter describing size of array for derivative orders *****/
    ndorder    = NDORDER;

    /**** Parameter describing interpolation results storage *****/
    rhint      = DF_MATRIX_STORAGE_ROWS;

    /***** Generate array of uniformly distributed break points *****/
    errcode = dUniformRandSortedData( x, left, right, nx );
    CheckDfError(errcode);

    /***** Generate function y = sin(2 * Pi * freq * x) *****/
    errcode = dSinDataNotUniformGrid( y, x, FREQ, nx );
    CheckDfError(errcode);

    /***** Generate interpolation sites *****/
    for (i = 0; i < N; i++)
    {
        site[i] = x[i];
    }

    /***** Create Data Fitting task *****/
    errcode = dfdNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Edit task parameters for look up interpolant *****/
    errcode = dfdEditPPSpline1D( task, sorder, stype, 0, 0, 0, 0, 0, 0 );
    CheckDfError(errcode);

    /***** Interpolate using lookup method *****/
    errcode = dfdInterpolate1D( task, DF_INTERP, DF_METHOD_PP,
                                    nsite, site, sitehint, ndorder,
                                    dorder, datahint, r, rhint, 0 );
    CheckDfError(errcode);

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Check interpolation results *****/
    for( j = 0; j < N; j++ )
    {
        if ( DF_ABS( r[j] - y[j] ) > EPSILON_DOUBLE )
            errnums++;
    }

    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);
    printf("Number of sites : %d\n", (int)nsite);

    /***** rint given function and computed results *****/
    printf("\n i  x(i)        y(i)       r(i)\n");
    for( j = 0; j < nx; j++ )
    {
        printf(" %+lf   %+lf   %+lf\n", x[j], y[j], r[j]);
    }

    /***** Print summary of the test *****/
    if (errnums != 0)
    {
        printf("\n\nError: Computed interpolation results");
        printf(" are incorrect\n");
        return 1;
    }
    else
    {
        printf("\n\nComputed interpolation results");
        printf(" are correct\n");
    }

    return 0;
}
