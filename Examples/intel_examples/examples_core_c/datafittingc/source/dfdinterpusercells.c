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
!    Construction and interpolation of linear spline using user provided
!  cell indices Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N              7 // number of break points
#define NY             1 // number of functions
#define NSITE         10 // number of interpolation sites

#define NSCOEFF      DF_PP_LINEAR*(N-1)  // number of spline coefficients

#define LLIM_X      -10.0 // left  limit of interpolation interval
#define RLIM_X       10.0 // right limit of interpolation interval

#define FREQ         1.75

int main()
{
    DFTaskPtr task;                     // Data Fitting task descriptor
    MKL_INT sorder;                     // spline order
    MKL_INT stype;                      // spline type
    MKL_INT nx;                         // number of break points
    MKL_INT ny;                         // number of functions
    MKL_INT xhint;                      // additional info about break points
    MKL_INT yhint;                      // additional info about function
    MKL_INT scoeffhint;                 // additional info about spline
                                        // coefficients
    MKL_INT bc_type;                    // boundary conditions type
    MKL_INT ic_type;                    // internal conditions type

    MKL_INT nsite;                      // number of interpolation sites
    MKL_INT sitehint;                   // additional info about interpolation
                                        // sites
    MKL_INT ndorder = 0;                // size of array describing derivative
                                        // orders
    MKL_INT *dorder = 0;                // spline values will be computed
    MKL_INT rhint;                      // interpolation results storage format

    double left  = LLIM_X;
    double right = RLIM_X;
    double x[N];                        // array of break points
    double y[N];                        // function values
    double *bc;
    double *ic;
    double scoeff[NSCOEFF];             // array of spline coefficients
    double site[NSITE];                 // interpolation sites
    double r[NSITE];                    // interpolation results
    double *datahint;                   // additional info about structure
                                        // of arrays x and y
    MKL_INT cell[NSITE];                // cell search results

    double ref_r[NSITE];                // reference interpolation results
    MKL_INT ref_cell[NSITE];            // reference cell search results

    int i, errcode, errnums;

    errcode = DF_STATUS_OK;
    errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/
    sorder = DF_PP_LINEAR;
    stype  = DF_PP_DEFAULT;

    /***** Parameters describing interpolation interval *****/
    nx          = N;
    xhint       = DF_NON_UNIFORM_PARTITION;

    /***** Parameters describing function *****/
    ny          = NY;
    yhint       = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    /* No boundary conditions are provided for linear spline */
    bc_type    = DF_NO_BC;
    bc         = 0;

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for linear spline */
    ic_type    = DF_NO_IC;
    ic         = 0;

    /***** Parameters describing spline coefficients storage *****/
    scoeffhint  = DF_NO_HINT;

    /***** Parameters describing interpolation sites *****/
    nsite      = NSITE;
    sitehint   = DF_NON_UNIFORM_PARTITION;

    /***** Additional info about structure of arrays x and y *****/
    /* No additional info is provided */
    datahint   = 0;

    /**** Parameter describing interpolation results storage *****/
    rhint      = DF_NO_HINT;

    /***** Generate array of uniformly distributed break points *****/
    errcode = dUniformRandSortedData( x, left, right, nx );
    CheckDfError(errcode);

    /***** Generate interpolation sites *****/
    errcode = dUniformRandData( site, left, right, nsite );
    CheckDfError(errcode);

    /***** Create Data Fitting task *****/
    errcode = dfdNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Perform cells search *****/
    errcode = dfdSearchCells1D( task, DF_METHOD_STD, nsite, site,
                                    sitehint, datahint, cell );
    CheckDfError(errcode);

    /***** Edit task parameters for linear spline construction *****/
    errcode = dfdEditPPSpline1D( task, sorder, stype, bc_type, bc,
                                 ic_type, ic, scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Generate functions y = sin(2 * Pi * freq * x) *****/
    errcode = dSinDataNotUniformGrid( y, x, FREQ, nx );
    CheckDfError(errcode);

    /***** Construct linear spline using STD method *****/
    errcode = dfdConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);

    /***** Interpolate using user provided cell indices *****/
    errcode = dfdInterpolate1D( task, DF_INTERP_USER_CELL, DF_METHOD_PP,
                                nsite, site, sitehint, ndorder,
                                dorder, datahint, r, rhint, cell );
    CheckDfError(errcode);

    /***** Check cell search results *****/
    errcode = dFindCells( nx, x, nsite, site, ref_cell );

    for ( i = 0; i < nsite; i++ )
    {
        if ( ref_cell[i] != cell[i] ) errnums++;
    }

    /***** Check interpolation results *****/
    errcode = dCheckLinInterpRes( nx, x, ny, scoeff, nsite, site,
                                  ndorder, dorder, r, ref_r );
    if (errcode < 0) errnums++;

    /***** Print results *****/
    printf("Number of break points:        %d\n", (int)nx);
    printf("Number of interpolation sites: %d\n", (int)nsite);

    /***** Print given function *****/
    printf(" i    x(i)       y(i)\n");
    for (i = 0; i < nx; i++)
    {
        printf(" %1d   %+.6lf  %+.6lf\n", i, x[i], y[i]);
    }

    /***** Print cell search and interpolation results *****/
    printf("\nResults of cell search and linear interpolation:\n");
    printf("  Site      Cell index      Interpolation result\n");
    printf("                            Obtained         Expected\n");
    for (i = 0; i < nsite; i++)
        printf(" %+.6lf    %8d     %+.6lf     %+.6lf\n", site[i], (int)cell[i],
            r[i], ref_r[i]);
    printf("\n");

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Print summary of the test *****/
    if (errnums != 0)
    {
        printf("\nError: Computed interpolation results");
        printf(" are incorrect\n");
        return 1;
    }
    else
    {
        printf("\nComputed interpolation results");
        printf(" are correct\n");
    }

    return 0;
}
