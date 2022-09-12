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
!    Calculation of Bessel cubic spline coefficients and spline evaluation
!   with user-defined extrapolation Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N              7 // number of breakpoints
#define NY             1 // number of datasets to interpolate
#define NSITE         10 // number of interpolation sites
#define NDORDER        1 // size of array describing derivative orders
                         // to compute

#define NSCOEFF        NY*(N-1)*DF_PP_CUBIC

#define LLIM_X         0.0  // left  limit of interpolation interval
#define RLIM_X         3.0  // right limit of interpolation interval
#define LLIM_SITE     -1.0  // left  limit of interpolation sites
#define RLIM_SITE      4.0  // right limit of interpolation sites
#define FREQ           0.7

/* Structure containing linear extrapolation parameters */
typedef struct _extrap_params
{
    double x;
    double y;
    double secant;
} extrap_params;

/*******************************************************************************
!   Definition of the call back for linear extrapolation.
!
! API
!   int  linear_extrap( MKL_INT64* n, MKL_INT64 cell[], double site[],
!                       double r[], void *user_params,
!                       dfInterpCallBackLibraryParams* library_params )
!
! INPUT PARAMETERS:
!  n      - number of interpolation sites
!  cell   - array of size n containing indices of the cells to which the
!           interpolation sites in array 'site' belong
!  site   - array of size n that holds the interpolation sites
!  user_params - pointer to the structure contatining extrapolation parameters
!  library_params - pointer to library-defined parameters
!
! OUTPUT PARAMETERS:
!  r      - array of integration results
!
! RETURN VALUE:
!  The status returned by the callback function:
!   0  - indicates successful completion of the callback operation
!   <0 - error
!   >0 - warning
*******************************************************************************/
int  linear_extrap( MKL_INT64* n, MKL_INT64 cell[], double site[], double r[],
                    void *user_params,
                    dfInterpCallBackLibraryParams* library_params );

int main()
{
    DFTaskPtr task;                       // Data Fitting task descriptor
    MKL_INT sorder;                       // spline order
    MKL_INT stype;                        // spline type
    MKL_INT nx;                           // number of break points
    MKL_INT xhint;                        // additional info about break points
    MKL_INT ny;                           // number of functions
    MKL_INT yhint;                        // functions storage format
    MKL_INT scoeffhint;                   // additional info about spline
                                          // coefficients
    MKL_INT bc_type;                      // boundary conditions type
    MKL_INT ic_type;                      // internal conditions type
    MKL_INT nsite;                        // number of interpolation sites
    MKL_INT sitehint;                     // additional info about interpolation
                                          // sites
    MKL_INT ndorder;                      // size of array describing derivative
                                          // orders
    MKL_INT dorder = 1;                   // only spline value will be computed
    MKL_INT rhint;                        // interpolation results storage format
    MKL_INT *cell_idx;                    // indices of cells containing
                                          // interpolation sites
    double left = LLIM_X, right = RLIM_X; // limits of interpolation interval
    double lsite = LLIM_SITE;             // left  limit of interpolation sites
    double rsite = RLIM_SITE;             // right limit of interpolation sites
    double x[N];                          // break points
    double y[N*NY];                       // function values
    double *bc;                           // boundary condition
    double *ic;                           // internal condition
    double scoeff[NSCOEFF];               // Bessel spline coefficients
    double site[NSITE];                   // array of interpolation sites
    double r[NDORDER*NSITE];              // spline evaluation results
    double *datahint;                     // additional info about structure
                                          // of arrays x and y
    double *y_cur, *scoeff_cur;

    dfdInterpCallBack le_cb, re_cb;       // interpolation call backs
    extrap_params le_params, re_params;   // interpolation call backs parameters

    double left_val[NY*(N-1)], right_val[NY*(N-1)];
    double left_der1[NY*(N-1)], right_der1[NY*(N-1)];
    double ref_r[NDORDER*NSITE];

    double freq;
    MKL_INT64 tmp_cell;
    MKL_INT64 n1 = 1;
    int i, j, errcode = 0;
    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/

    /***** Parameters describing order and type of the spline *****/
    sorder     = DF_PP_CUBIC;
    stype      = DF_PP_BESSEL;

    /***** Parameters describing interpolation interval *****/
    nx         = N;
    xhint      = DF_NON_UNIFORM_PARTITION;

    /***** Parameters describing functions *****/
    ny         = NY;
    yhint      = DF_MATRIX_STORAGE_ROWS;

    /***** Parameter describing additional info about spline
           coefficients *****/
    scoeffhint   = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    /* No additional parameters are provided for Not-a-Knot
       boundary conditions */
    bc_type    = DF_BC_NOT_A_KNOT;
    bc = 0;

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for Bessel cubic spline */
    ic_type    = DF_NO_IC;
    ic         = 0;

    /***** Parameters describing interpolation sites *****/
    nsite      = NSITE;
    sitehint   = DF_SORTED_DATA;

    /***** Additional info about structure of arrays x and y *****/
    /* No additional info is provided */
    datahint   = 0;

    /**** Parameter describing interpolation results storage *****/
    rhint      = DF_NO_HINT;

    /**** Parameter describing array of cell indices *****/
    /* No cell indices are computed */
    cell_idx   = 0;

    /**** Parameter describing size of array for derivative orders *****/
    ndorder = NDORDER;

    /***** Generate independent variables array with
           quasi-uniform break points *****/
    errcode = dQuasiUniformData( x, left, right, nx );
    CheckDfError(errcode);

    /***** Generate interpolation sites *****/
    errcode = dUniformRandSortedData( site, lsite, rsite, nsite );
    CheckDfError(errcode);

    /***** Generate function y = sin(2 * Pi * freq * x) *****/
    freq = FREQ;
    errcode = dSinDataNotUniformGrid( y, x, freq, nx );
    CheckDfError(errcode);

    /***** Call backs for integration on the outer intervals *****/
    le_cb = linear_extrap;
    re_cb = linear_extrap;

    /***** Set left call back parameters ******/
    le_params.x = x[0];
    le_params.y = y[0];
    le_params.secant = (y[1] - y[0]) / (x[1] - x[0]);

    /***** Set right call back parameters ******/
    re_params.x = x[nx-2];
    re_params.y = y[nx-2];
    re_params.secant = (y[nx-1] - y[nx-2]) / (x[nx-1] - x[nx-2]);

    /***** Create Data Fitting task *****/
    errcode = dfdNewTask1D( &task, nx, x, xhint, ny, y, DF_NO_HINT );
    CheckDfError(errcode);

    /***** Edit task parameters for Bessel cubic spline coefficients
           computation *****/
    errcode = dfdEditPPSpline1D( task, sorder, stype, bc_type, 0, ic_type, 0,
                                 scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Compute Bessel cubic spline coefficients using STD method *****/
    errcode = dfdConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);

    /***** Interpolate and use call backs for the left and right
           extrapolation *****/
    errcode = dfdInterpolateEx1D( task, DF_INTERP, DF_METHOD_PP, nsite, site,
                                  sitehint, ndorder, &dorder, datahint,
                                  r, rhint, cell_idx, le_cb, &le_params,
                                  re_cb, &re_params, 0, 0, 0, 0 );
    CheckDfError(errcode);

    /***** Check computed coefficients *****/

    /***** Check spline values in break points *****/
    errcode = dCheckCubBreakPoints( nx, x, ny, y, scoeff, left_val, right_val );
    if ( errcode < 0 ) errnums++;

    /***** Check that spline 1st derivatives are equal for left
           and right piece of the spline for each break point *****/
    errcode = dCheckCub1stDerConsistency( nx, x, ny, scoeff,
                                          left_der1, right_der1 );
    if ( errcode < 0 ) errnums++;

    /***** Check computed interpolation results *****/
    errcode = dCheckCubInterpRes( nx, x, ny, scoeff,
                                  nsite, site, ndorder, &dorder,
                                  r, ref_r );
    if ( errcode < 0 ) errnums++;

    /***** Check computed extrapolation results *****/

    for ( i = 0; i < nsite; i++ )
    {
        if ( site[i] < x[0] )
        {
            /* Check left extrapolation results */

            dFindCells( nx, x, (MKL_INT)1, &site[i], (MKL_INT*)(&tmp_cell) );

            errcode = linear_extrap( &n1, &tmp_cell, &site[i], &ref_r[i],
                &le_params, NULL );

            if ( DF_ABS(ref_r[i] - r[i]) > EPSILON_DOUBLE )
                errnums++;
        }
        if ( site[i] > x[nx-1] )
        {
            /* Check right extrapolation results */

            dFindCells( nx, x, (MKL_INT)1, &site[i], (MKL_INT*)(&tmp_cell) );
            tmp_cell--;

            errcode = linear_extrap( &n1, &tmp_cell, &site[i], &ref_r[i],
                &re_params, NULL );

            if ( DF_ABS(ref_r[i] - r[i]) > EPSILON_DOUBLE )
                errnums++;
        }
    }

    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);

    /***** Print given tabular function *****/
    printf("\n i  x(i)        y(i)\n");

    for( j = 0; j < nx; j++ )
    {
        printf(" %1d %+lf   %+lf\n", j, x[j], y[j]);
    }

    /***** Print computed spline coefficients *****/

    printf("\nCoefficients are calculated for a polynomial of the form:\n\n");
    printf("Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 + Di*(x - x(i))^3\n");
    printf("    where x(i) <= x < x(i+1)\n");
    printf("\nSpline coefficients for Y:\n");
    printf(" i    Ai            Bi            Ci            Di      ");
    printf("      Pi(x(i))      Pi(x(i+1))    Pi'(x(i))     Pi'(x(i+1))\n");

    for( j = 0; j < nx-1; j++ )
    {
        printf(" %1d %+11.6f   %+11.6f   %+11.6f   %+11.6f   %+11.6f   %+11.6f",
                j, scoeff[j*sorder], scoeff[j*sorder + 1],
                scoeff[j*sorder + 2], scoeff[j*sorder + 3],
                right_val[j], left_val[j]);
        printf("   %+11.6f   %+11.6f\n", right_der1[j], left_der1[j]);
    }

    /***** Print interpolation results ******/
    printf("\nResults of interpolation:\n");
    printf("    Sites         Obtained     Expected\n");
    for ( i = 0; i < nsite; i++ )
    {
        printf(" %+11.6lf  %11.6lf  %11.6lf\n", site[i], r[i], ref_r[i]);
    }

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Printing summary of the test *****/
    if (errnums != 0) {
        printf("\n\nError: Computed Bessel cubic spline coefficients");
        printf(" or spline values are incorrect\n");
        return 1;
    }
    else {
        printf("\n\nComputed Bessel cubic spline coefficients");
        printf(" and spline values are correct\n");
    }

    return 0;
}

/*******************************************************************************
!   Call back for linear extrapolation.
!
! API
!   int  linear_extrap( MKL_INT* n, MKL_INT cell[], double site[],
!                       double r[], void *user_params,
!                       dfInterpCallBackLibraryParams* library_params )
!
! INPUT PARAMETERS:
!  n      - number of interpolation sites
!  cell   - array of size n containing indices of the cells to which the
!           interpolation sites in array 'site' belong
!  site   - array of size n that holds the interpolation sites
!  user_params - pointer to the structure contatining extrapolation parameters
!  library_params - pointer to library-defined parameters
!
! OUTPUT PARAMETERS:
!  r      - array of integration results
!
! RETURN VALUE:
!  The status returned by the callback function:
!   0  - indicates successful completion of the callback operation
!   <0 - error
!   >0 - warning
*******************************************************************************/
int  linear_extrap( MKL_INT64* n, MKL_INT64 cell[], double site[], double r[],
                    void *user_params,
                    dfInterpCallBackLibraryParams* library_params )
{
    MKL_INT64 i;
    extrap_params *p = (extrap_params*)user_params;

    for ( i = 0; i < n[0]; i++ )
    {
        r[i] = p->y + p->secant * ( site[i] - p->x );
    }

    return 0;
}
