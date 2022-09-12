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
!    Construction of Hermite cubic spline and interpolation
!  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N              7 // number of break points
#define NY             1 // number of functions
#define NBC            2 // number of boundary conditions
#define NIC          N-2 // number of internal conditions
#define NSITE         15 // number of interpolation sites
#define NDORDER        1 // size of array describing derivative orders
                         // to compute

#define NSCOEFF     (NY*(N-1)*DF_PP_CUBIC)   // total number of spline
                                             // coefficients

#define LLIM_X       0.0 // left  limit of interpolation interval
#define RLIM_X       3.0 // right limit of interpolation interval
#define FREQ         0.3

int main()
{
    DFTaskPtr task;                     // Data Fitting task descriptor
    MKL_INT sorder;                     // spline order
    MKL_INT stype;                      // spline type
    MKL_INT nx;                         // number of break points
    MKL_INT xhint;                      // additional info about break points
    MKL_INT ny;                         // number of functions
    MKL_INT yhint;                      // additional info about function
    MKL_INT nscoeff;                    // number of spline coefficients
    MKL_INT scoeffhint;                 // additional info about spline
                                        // coefficients
    MKL_INT bc_type;                    // boundary conditions type
    MKL_INT nbc;                        // number of boundary conditions
    MKL_INT ic_type;                    // internal conditions type
    MKL_INT nic;                        // number of internal conditions
    MKL_INT nsite;                      // number of interpolation sites
    MKL_INT sitehint;                   // additional info about interpolation
                                        // sites
    MKL_INT ndorder;                    // size of array describing derivative
                                        // orders
    MKL_INT dorder[] = {1};             // spline values will be computed
    MKL_INT rhint;                      // interpolation results storage format
    MKL_INT *cell_idx;                  // indices of cells containing
                                        // interpolation sites
    float left  = LLIM_X;               // left limit of the interpolation
                                        // interval
    float right = RLIM_X;               // right limit of the interpolation
                                        // interval
    float x[2];                         // limits of the interpolation interval
    float xx[N];                        // partition in full format
    float y[NY*N];                      // function values
    float ic[NIC];                      // array of internal conditions
    float bc[NBC];                      // array of boundary conditions
    float scoeff[NSCOEFF];              // array of spline coefficients
    float site[NSITE];                  // array of interpolation sites
    float r[NDORDER*NSITE];             // spline evaluation results
    float *datahint;                    // additional info about structure
                                        // of arrays x and y

    float freq, delta;
    float left_val[N-1], right_val[N-1];
    float left_der[N-1], right_der[N-1];
    float ref_r[(NDORDER+1)*NSITE];

    int i, j, errcode = 0;
    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/

    sorder = DF_PP_CUBIC;
    stype  = DF_PP_HERMITE;

    /***** Parameters describing interpolation interval *****/
    nx          = N;
    xhint       = DF_UNIFORM_PARTITION;
    /* Limits of interpolation interval are provided in case
       of uniform partition */
    x[0]        = left;
    x[1]        = right;

    /***** Parameters describing function *****/
    ny          = NY;
    yhint       = DF_NO_HINT;

    /***** Parameters describing spline coefficients storage *****/
    nscoeff     = NSCOEFF;
    scoeffhint  = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    bc_type     = DF_BC_1ST_LEFT_DER | DF_BC_1ST_RIGHT_DER;

    /***** Parameters describing internal conditions *****/
    ic_type     = DF_IC_1ST_DER;
    nic         = NIC;

    /***** Parameters describing interpolation sites *****/
    nsite      = NSITE;
    sitehint   = DF_NON_UNIFORM_PARTITION;

    /***** Additional info about structure of arrays x and y *****/
    /* No additional info is provided */
    datahint   = 0;

    /**** Parameter describing interpolation results storage *****/
    rhint      = DF_MATRIX_STORAGE_COLS;

    /**** Parameter describing array of cell indices *****/
    /* No cell indices are computed */
    cell_idx   = 0;

    /**** Parameter describing size of array for derivative orders *****/
    ndorder    = NDORDER;

    /***** Generate partition in full format *****/
    errcode = sUniformData( xx, left, right, nx );
    CheckDfError(errcode);

    /***** Generate function y = sin(2 * Pi * freq * x) *****/
    freq = FREQ;
    errcode = sSinDataNotUniformGrid( y, xx, freq, nx );
    CheckDfError(errcode);

    /***** Generate array of 1st derivatives needed for Hermite spline
           construction *****/
    errcode = sSinDerDataNotUniformGrid( ic, &xx[1], freq, nic );
    CheckDfError(errcode);

    /***** Generate boundary conditions *****/
    errcode = sSinDerDataNotUniformGrid( &bc[0], &xx[0],    freq, 1 );
    CheckDfError(errcode);
    errcode = sSinDerDataNotUniformGrid( &bc[1], &xx[nx-1], freq, 1 );
    CheckDfError(errcode);

    /***** Generate interpolation sites *****/
    errcode = sUniformRandData( site, left, right, nsite );
    CheckDfError(errcode);

    /***** Create Data Fitting task *****/
    errcode = dfsNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Edit task parameters for Hermite cubic spline construction *****/
    errcode = dfsEditPPSpline1D( task, sorder, stype, bc_type, bc,
                                 ic_type, ic, scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Construct Hermite cubic spline using STD method *****/
    errcode = dfsConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);

    /***** Interpolate using PP method *****/
    errcode = dfsInterpolate1D( task, DF_INTERP, DF_METHOD_PP,
                                nsite, site, sitehint, ndorder,
                                dorder, datahint, r, rhint, cell_idx );
    CheckDfError(errcode);

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Check computed coefficients *****/

    /***** Check spline values in break points *****/
    errcode = sCheckCubBreakPoints( nx, xx, ny, y, scoeff,
                                    left_val, right_val );
    if ( errcode < 0 ) errnums++;

    /***** Check that spline 1st derivatives are equal for left
           and right piece of the spline for each break point *****/
    errcode = sCheckCub1stDerConsistency( nx, xx, ny, scoeff,
                                          left_der, right_der );
    if ( errcode < 0 ) errnums++;

    /***** Check internal conditions *****/
    for( j = 0; j < nic; j++ )
    {
        if ( DF_ABS( ic[j] - left_der[j] ) > EPSILON_SINGLE )
            errnums++;
    }

    /***** Check boundary conditions *****/
    errcode = sCheckCubBC( nx, xx, ny, scoeff, bc_type, bc );
    if ( errcode < 0 ) errnums++;

    /***** Check results of interpolation *****/
    errcode = sCheckCubInterpRes( nx, xx, ny, scoeff, nsite, site, ndorder,
                                  dorder, r, ref_r );
    if (errcode < 0) errnums++;

    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);

    /***** Print given function *****/
    printf("\n i  x(i)        y(i)        y'(i)\n");

    printf(" %1d %+lf   %+lf   %+lf\n", 0, xx[0], y[0], bc[0]);
    for( j = 1; j < nx-1; j++ )
    {
        printf(" %1d %+lf   %+lf   %+lf\n", j,  xx[j], y[j], ic[j-1]);
    }
    printf(" %1d %+lf   %+lf   %+lf\n", (int)(nx-1),  xx[nx-1], y[nx-1], bc[1]);

    /***** Print computed spline coefficients *****/
    printf("\nCoefficients are calculated for a polynomial of the form:\n\n");
    printf("Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 + Di*(x - x(i))^3\n");
    printf("    where x(i) <= x < x(i+1)\n");
    printf("\nSpline coefficients for Y:\n");
    printf(" i    Ai            Bi            Ci            Di        ");
    printf("    P(x(i))       P(x(i+1)) ");
    printf("    P'(x(i))      P'(x(i+1))\n");

    for( j = 0; j < nx-1; j++ )
    {
        printf(" %1d %+11.6f   %+11.6f   %+11.6f   %+11.6f   %+11.6f   %+11.6f",
                j, scoeff[sorder*j], scoeff[sorder*j + 1], scoeff[sorder*j + 2],
                scoeff[sorder*j + 3], right_val[j], left_val[j]);
        printf("   %+11.6f   %+11.6f\n", right_der[j], left_der[j]);
    }

    printf("\n  i    Sites       Spline value\n");
    printf("                   Computed    Expected\n");
    for ( j = 0; j < nsite; j++ )
    {
        printf(" %2d %+11.6lf %+11.6lf %+11.6lf\n", j,
            site[j], r[j], ref_r[j]);
    }

    /***** Print summary of the test *****/
    if (errnums != 0)
    {
        printf("\n\nError: Computed Hermite cubic spline coefficients, ");
        printf(" spline values are incorrect\n");
        return 1;
    }
    else
    {
        printf("\n\nComputed Hermite cubic spline coefficients, ");
        printf(" spline values are correct\n");
    }

    return 0;
}
