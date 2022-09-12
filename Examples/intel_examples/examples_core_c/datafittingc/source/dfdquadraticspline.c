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
!    Construction and evaluation of quadratic spline Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N              9 // number of break points
#define NY             1 // number of functions
#define NSITE         10 // number of interpoaltion sites
#define NDORDER        1 // size of array describing derivative orders
                         // to compute

#define NSCOEFF        NY*(N-1)*DF_PP_QUADRATIC // total number of spline
                                                // coefficients

#define LLIM_X     -2.0 // left  limit of interpolation interval
#define RLIM_X      2.0 // right limit of interpolation interval
#define LLIM_SITE  -3.0 // left  limit of interpolation sites
#define RLIM_SITE   3.0 // right limit of interpolation sites
#define FREQ        0.5

int main()
{
    DFTaskPtr task;                     // Data Fitting task descriptor
    MKL_INT sorder;                     // spline order
    MKL_INT stype;                      // spline type
    MKL_INT nx;                         // number of break points
    MKL_INT xhint;                      // additional info about break points
    MKL_INT ny;                         // number of functions
    MKL_INT yhint;                      // additional info about function
    MKL_INT scoeffhint;                 // additional info about spline
                                        // coefficients
    MKL_INT bc_type;                    // boundary conditions type
    MKL_INT ic_type;                    // internal conditions type
    MKL_INT nsite;                      // number of interpolation sites
    MKL_INT sitehint;                   // additional info about interpolation
                                        // sites
    MKL_INT ndorder;                    // size of array describing derivative
                                        // orders
    MKL_INT dorder[] = {1};             // spline values will be computed
    MKL_INT rhint;                      // interpolation results storage format
    MKL_INT *cell_idx;                  // indices of cells containing
                                        // interpolation sites
    double x[N];                        // array of break points
    double y[N*NY];                     // function values
    double bc;                          // boundary condition
    double *ic;                         // internal conditions
    double scoeff[NSCOEFF];             // array of spline coefficients
    double site[2];                     // limits of interpolation sites
                                        // are provided if the sites are uniform
    double site_full[NSITE];            // array of interpolation sites
                                        // in full format
    double r[NDORDER*NSITE];            // spline evaluation results
    double *datahint;                   // additional info about structure
                                        // of arrays x and y

    double left = LLIM_X, right = RLIM_X;
    double freq;
    double left_val[N-1], right_val[N-1];
    double left_der[N-1], right_der[N-1];
    double ref_r[NDORDER*NSITE];

    int i, j, errcode = 0;
    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/

    sorder = DF_PP_QUADRATIC;
    stype  = DF_PP_DEFAULT;

    /***** Parameters describing interpolation interval *****/
    nx         = N;
    xhint      = DF_NON_UNIFORM_PARTITION;

    /***** Parameters describing function *****/
    ny         = NY;
    yhint      = DF_NO_HINT;

    /***** Parameters describing spline coefficients storage *****/
    scoeffhint = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    bc_type    = DF_BC_Q_VAL;
    bc         = 1.0;

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for quadratic spline */
    ic_type    = DF_NO_IC;
    ic         = 0;

    /***** Parameters describing interpolation sites *****/
    nsite      = NSITE;
    sitehint   = DF_UNIFORM_PARTITION;
    /* Limits of interpolation interval are provided if the sites are uniform */
    site[0]    = LLIM_SITE;
    site[1]    = RLIM_SITE;

    /***** Additional info about structure of arrays x and y *****/
    /* No additional info is provided */
    datahint   = 0;

    /**** Parameter describing interpolation results storage *****/
    rhint      = DF_NO_HINT;

    /**** Parameter describing array of cell indices *****/
    /* No cell indices are computed */
    cell_idx   = 0;

    /**** Parameter describing size of array for derivative orders *****/
    ndorder    = NDORDER;

    /***** Generate uniformly distributed break points *****/
    errcode = dUniformRandSortedData( x, left, right, nx );
    CheckDfError(errcode);

    /***** Generate function y = sin(2 * Pi * freq * x) *****/
    freq = FREQ;
    errcode = dSinDataNotUniformGrid( y, x, freq, nx );
    CheckDfError(errcode);

    /***** Create Data Fitting task *****/
    errcode = dfdNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Edit task parameters for quadratic spline construction *****/
    errcode = dfdEditPPSpline1D( task, sorder, stype, bc_type, &bc,
                                 ic_type, ic, scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Construct quadratic spline using STD method *****/
    errcode = dfdConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);

    /***** Interpolate using PP method *****/
    errcode = dfdInterpolate1D( task, DF_INTERP, DF_METHOD_PP, nsite, site,
                                sitehint, ndorder, dorder, datahint, r, rhint,
                                cell_idx );
    CheckDfError(errcode);

    /***** Check computed coefficients *****/

    /***** Check spline values in break points *****/
    errcode = dCheckQuadBreakPoints( nx, x, ny, y, scoeff,
                                     left_val, right_val );
    if ( errcode < 0 ) errnums++;

    /***** Check 1st derivatives in break points *****/
    errcode = dCheckQuad1stDerConsistency( nx, x, ny, scoeff,
                                           left_der, right_der );
    if ( errcode < 0 ) errnums++;

    /***** Check boundary conditions *****/
    errcode = dCheckQuadBC( nx, x, ny, scoeff, bc_type, &bc );
    if ( errcode < 0 ) errnums++;

    /***** Check results of interpolation *****/
    errcode = dUniformData( site_full, LLIM_SITE, RLIM_SITE, nsite );
    CheckDfError(errcode);

    errcode = dCheckQuadInterpRes( nx, x, ny, scoeff,
                                   nsite, site_full, ndorder, dorder,
                                   r, ref_r );
    if (errcode < 0) errnums++;

    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);

    /***** Print given function *****/
    printf("\n i  x(i)        y(i)\n");

    for( j = 0; j < nx; j++ )
    {
        printf(" %1d %+lf   %+lf\n", j, x[j], y[j]);
    }

    /***** Print computed spline coefficients *****/
    printf("\nCoefficients are calculated for a polynomial of the form:\n\n");
    printf("Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2\n");
    printf("    where x(i) <= x < x(i+1)\n");
    printf("\nSpline coefficients for Y:\n");
    printf(" i    Ai            Bi            Ci        ");
    printf("    P(x(i))       P(x(i+1)) ");
    printf("    P'(x(i))      P'(x(i+1))\n");

    for( j = 0; j < nx-1; j++ )
    {
        printf(" %1d %+11.6f   %+11.6f   %+11.6f   %+11.6f   %+11.6f", j,
                scoeff[j*sorder], scoeff[j*sorder + 1], scoeff[j*sorder + 2],
                right_val[j], left_val[j]);
        printf("   %+11.6f   %+11.6f\n", right_der[j], left_der[j]);
    }

    /***** Print interpolation results ******/
    printf("\nResults of interpolation:\n");
    printf("                    Function value\n");
    printf("    Sites         Obtained     Expected\n");
    for ( i = 0; i < nsite; i++ )
    {
        printf(" %+11.6lf", site_full[i]);
        printf("   %+11.6lf  %+11.6lf\n", r[i],
                ref_r[i]);
    }

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Print summary of the test *****/
    if (errnums != 0)
    {
        printf("\n\nError: Computed quadratic spline coefficients");
        printf(" or values are incorrect\n");
        return 1;
    }
    else
    {
        printf("\n\nComputed quadratic spline coefficients and values ");
        printf(" are correct\n");
    }

    return 0;
}
