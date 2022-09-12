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
!    Construction and integration of natural cubic spline
!  Example Program Text
!******************************************************************************/

#include <stdio.h>

#include "mkl.h"
#include "errcheck.inc"
#include "generatedata.inc"
#include "rescheck.inc"

#define N              6 // number of break points
#define NY             1 // number of functions
#define NLIM           4 // number of pairs of integration limits

#define NSCOEFF        (NY*(N-1)*DF_PP_CUBIC)  // total number of spline
                                               // coefficients

#define LLIM_X       0.0 // left limit of interpolation interval
#define RLIM_X       7.0 // right limit of interpolation interval
#define LEFT_LLIM    1.0
#define RIGHT_LLIM   3.0
#define LEFT_RLIM    4.0
#define RIGHT_RLIM   5.5
#define FREQ        0.15


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
    MKL_INT nlim;                       // number of pairs of integration limits
    MKL_INT llimhint;                   // additional info about the structure
                                        // of left integration limits
    MKL_INT rlimhint;                   // additional info about the structure
                                        // of right integration limits
    MKL_INT rhint;                      // integration results storage format
    double x[2];                        // limits of the interpolation interval
    double y[NY*N];                     // function values
    double *ic;                         // internal conditions
    double *bc;                         // boundary conditions
    double scoeff[NSCOEFF];             // array of spline coefficients
    double llim[NLIM];                  // left  integration limits
    double rlim[NLIM];                  // right integration limits
    double *ldatahint, *rdatahint;      // additional info about the structure
                                        // of x and integration limits
    double r[NY*NLIM];                  // integration results
    double ref_r[NY*NLIM];              // reference integration results

    double left = LLIM_X, right = RLIM_X;
    double freq;
    double xx[N];
    double left_val[N-1], right_val[N-1];
    double left_der1[N-1], right_der1[N-1];
    double left_der2[N-1], right_der2[N-1];

    int i, j, errcode = 0;
    int errnums = 0;

    /***** Initializing parameters for Data Fitting task *****/

    sorder     = DF_PP_CUBIC;
    stype      = DF_PP_NATURAL;

    /***** Parameters describing interpolation interval *****/
    nx         = N;
    xhint      = DF_UNIFORM_PARTITION;
    /* Limits of interpolation interval are provided in case of uniform grid */
    x[0]       = left;
    x[1]       = right;

    /***** Parameters describing function *****/
    ny         = NY;
    yhint      = DF_NO_HINT;

    /***** Parameters describing spline coefficients storage *****/
    nscoeff    = NSCOEFF;
    scoeffhint = DF_NO_HINT;

    /***** Parameters describing boundary conditions type *****/
    bc_type    = DF_BC_FREE_END;
    /* No additional parameters are provided for Free-End boundary conditions */
    bc         = 0;

    /***** Parameters describing internal conditions type *****/
    /* No internal conditions are provided for natural cubic spline */
    ic_type    = DF_NO_IC;
    ic         = 0;

    /***** Parameters decsribing integration limits *****/
    nlim = NLIM;
    llimhint = DF_NON_UNIFORM_PARTITION;
    rlimhint = DF_NON_UNIFORM_PARTITION;

    /***** Additional information about the structure of integration
           limits *****/
    /* No additional info is provided */
    ldatahint = 0;
    rdatahint = 0;

    /***** Parameter dascribing integration results storage format *****/
    rhint = DF_NO_HINT;

    /***** Generate function y = sin(2 * Pi * freq * x) *****/
    errcode = dSinDataUniformGrid( y, FREQ, left, right, nx );
    CheckDfError(errcode);

    /***** Generate limits of integration intervals *****/
    errcode = dUniformRandSortedData( llim, LEFT_LLIM, RIGHT_LLIM, nlim );
    CheckDfError(errcode);
    errcode = dUniformRandSortedData( rlim, LEFT_RLIM, RIGHT_RLIM, nlim );
    CheckDfError(errcode);

    /***** Create Data Fitting task *****/
    errcode = dfdNewTask1D( &task, nx, x, xhint, ny, y, yhint );
    CheckDfError(errcode);

    /***** Edit task parameters for natural cubic spline construction *****/
    errcode = dfdEditPPSpline1D( task, sorder, stype, bc_type, bc,
                                 ic_type, ic, scoeff, scoeffhint );
    CheckDfError(errcode);

    /***** Construct natural cubic spline using STD method *****/
    errcode =  dfdConstruct1D( task, DF_PP_SPLINE, DF_METHOD_STD );
    CheckDfError(errcode);

    /***** Compute integrals *****/
    errcode = dfdIntegrate1D( task, DF_METHOD_PP, nlim, llim, llimhint,
                              rlim, rlimhint, ldatahint, rdatahint,
                              r, rhint );
    CheckDfError(errcode);

    /***** Check computed coefficients *****/
    errcode = dUniformData( xx, left, right, nx );
    CheckDfError(errcode);

    /***** Check spline values in break points *****/
    errcode = dCheckCubBreakPoints( nx, xx, ny, y, scoeff,
                                    left_val, right_val );
    if ( errcode < 0 ) errnums++;

    /***** Check that spline 1st derivatives are equal for left
           and right piece of the spline for each break point *****/
    errcode = dCheckCub1stDerConsistency( nx, xx, ny, scoeff,
                                          left_der1, right_der1 );
    if ( errcode < 0 ) errnums++;

    /***** Check that spline 2nd derivatives are equal for left
           and right piece of the spline for each break point *****/
    errcode = dCheckCub2ndDerConsistency( nx, xx, ny, scoeff,
                                          left_der2, right_der2 );
    if ( errcode < 0 ) errnums++;

    /***** Check boundary conditions *****/
    errcode = dCheckCubBC( nx, xx, ny, scoeff, bc_type, bc );
    if ( errcode < 0 ) errnums++;

    /***** Check integration results *****/
    errcode = dCheckCubIntegrRes( nx, xx, ny, scoeff, nlim, llim, rlim, r, ref_r );
    if ( errcode < 0 ) errnums++;

    /***** Print results *****/
    printf("Number of break points : %d\n", (int)nx);

    /***** Print given function *****/
    printf("\n  X           Y\n");

    for( j = 0; j < nx; j++ )
    {
        printf(" %+lf   %+lf\n", xx[j], y[j]);
    }

    /***** Print computed spline coefficients *****/
    printf("\nCoefficients are calculated for a polynomial of the form:\n\n");
    printf("Pi(x) = Ai + Bi*(x - x(i)) + Ci*(x - x(i))^2 + Di*(x - x(i))^3\n");
    printf("    where x(i) <= x < x(i+1)\n");
    printf("\nSpline coefficients for Y:\n");
    printf(" i    Ai            Bi            Ci            Di        ");
    printf("    P(X[i])       P(X[i+1]) ");
    printf("    P'(X[i])      P'(X[i+1])    P\"(X[i])      P\"(X[i+1])\n");

    for( j = 0; j < nx-1; j++ )
    {
        printf(" %1d %+11.6f   %+11.6f   %+11.6f   %+11.6f   %+11.6f   %+11.6f",
                j, scoeff[j*sorder], scoeff[j*sorder + 1],
                scoeff[j*sorder + 2], scoeff[j*sorder + 3],
                right_val[j], left_val[j]);
        printf("   %+11.6f   %+11.6f   %+11.6f   %+11.6f\n",
                right_der1[j], left_der1[j], right_der2[j], left_der2[j]);
    }

    printf("\nIntegration results for Y:\n");
    printf("Integration interval      Result\n");
    for ( j = 0; j < nlim; j++ )
    {
        printf(" ( %+4.1lf, %+4.1lf )        %+11.6lf\n", llim[j], rlim[j],
            r[j]);
    }

    /***** Delete Data Fitting task *****/
    errcode = dfDeleteTask( &task );
    CheckDfError(errcode);

    /***** Print summary of the test *****/
    if (errnums != 0)
    {
        printf("\n\nError: Computed natural cubic spline coefficients");
        printf(" or integrals are incorrect\n");
        return 1;
    }
    else
    {
        printf("\n\nComputed natural cubic spline coefficients");
        printf(" and integrals are correct\n");
    }

    return 0;
}
